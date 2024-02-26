/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_MULTIPLICATIVEINFLATION_H_
#define OOPS_BASE_MULTIPLICATIVEINFLATION_H_

#include <functional>
#include <memory>
#include <string>
#include <vector>

#include "atlas/field.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/InflationBase.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/State.h"
#include "oops/base/StateSet.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

namespace oops {

class MultiplicativeParameters : public InflationParameters{
  OOPS_CONCRETE_PARAMETERS(MultiplicativeParameters, InflationParameters);
 public:
    OptionalParameter<double> factor{"factor", this};
    OptionalParameter<std::vector<double>> levels{"levels", this};
    OptionalParameter<std::vector<double>> factors{"factors", this};
};

/// \brief Class for scaling the analysis perturbation by a factor
///
/// \details This class can be used for either analysis increments or states.
///
/// For an ensemble of states, the perturbation from the analysis mean
/// is multiplied by the factor specified in the configuration
///
/// For an ensemble of analysis increments:
///
/// \f$ dxa_i <- (\alpha * dxa_i) + (1 - \alpha) * (xb'_i - dxa_mean)\f$,
///
/// where dxa_i refers to the analysis increment for member i,
/// dxa_mean is the mean field of the ensemble of increments,
/// xb'_i is the background perturbation from the mean
/// and \alpha is the prescribed factor.
///
/// The factor can be applied to the entire analysis state, or to specific model levels
/// For model levels, a vector of factors and a vector of levels must be specified.
/// At each of the levels specified, the corresponding factor is applied here and then
/// interpolated between this and the neighbouring levels given.
/// If the top and bottom level are not defined in the input vector, they are set to
/// have a factor of 1.0
///
/// See section 3:
/// Inverarity, G.W., Tennant, W.J., Anton, L., Bowler, N.E., Clayton, A.M., Jardak, M.,
/// et al. (2023) Met Office MOGREPS-G initialisation using an ensemble of hybrid
/// four-dimensional ensemble variational (En-4DEnVar) data assimilations.
/// Quarterly Journal of the Royal Meteorological Society, 149(753), 1138–1164.
/// Available from: https://doi.org/10.1002/qj.4431

template <typename MODEL> class MultiplicativeInflation : public InflationBase<MODEL> {
  typedef Geometry<MODEL>                   Geometry_;
  typedef Increment<MODEL>                  Increment_;
  typedef IncrementSet<MODEL>               IncrementSet_;
  typedef State<MODEL>                      State_;
  typedef StateSet<MODEL>                   StateSet_;

 public:
  typedef MultiplicativeParameters          Parameters_;
  MultiplicativeInflation(const Parameters_ &, const Geometry_ &,
                          const StateSet_ &, const Variables &);

  void doInflation(IncrementSet_ &) override;
  void doInflation(StateSet_ &) override;

 private:
  void inflateByLevel(IncrementSet_ &);
  void inflateByLevel(StateSet_ &);
  void interpolateFactors(std::vector<double> &, atlas::FieldSet &);
  const Parameters_ params_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
MultiplicativeInflation<MODEL>::MultiplicativeInflation(const Parameters_ & params,
                                                        const Geometry_ & geom,
                                                        const StateSet_ & bens,
                                                        const Variables & vars)
  : InflationBase<MODEL>(geom, bens, vars), params_(params)
{
  Log::trace() << "MultiplicativeInflation::set up MultiplicativeInflation" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void MultiplicativeInflation<MODEL>::doInflation(IncrementSet_ & anEns) {
    Log::test() << "Multiplicative Background member 1:" << this->background()[0] << std::endl;
    Log::test() << "Multiplicative Analysis Increment member 1:" << anEns[0] << std::endl;

  if (params_.levels.value() != boost::none && params_.factors.value() != boost::none &&
      params_.factor.value() == boost::none) {
    Log::trace() << "MultiplicativeInflation::Increment inflate by level starting" << std::endl;
    this->inflateByLevel(anEns);
    Log::trace() << "MultiplicativeInflation::Increment inflation by level done" << std::endl;
  } else if (params_.factor.value() != boost::none && params_.levels.value() == boost::none &&
             params_.factors.value() == boost::none) {
    Log::trace() << "MultiplicativeInflation::Increment inflate all levels starting" << std::endl;
    const double & factor =  *params_.factor.value();
    StateSet_ bgMean = this->background().ens_mean();
    IncrementSet_ anMean = anEns.ens_mean();

    Increment_ bgPert(this->geometry(), this->vars(), this->background()[0].validTime());

    for (size_t ii = 0; ii < anEns.size(); ++ii) {
      anEns[ii] *= factor;

      bgPert.diff(this->background()[ii], bgMean[0]);

      bgPert -= anMean[0];
      bgPert *= (factor - 1.0);
      anEns[ii] += bgPert;
    }
    Log::trace() << "MultiplicativeInflation::Increment inflate all levels done" << std::endl;
  } else {
    ABORT("Incorrect input settings");
  }
  Log::test() << "Multiplicative Updated Analysis Increment member 1:" << anEns[0] << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void MultiplicativeInflation<MODEL>::doInflation(StateSet_ & anEns) {
  Log::test() << "Multiplicative Background member 1:" << this->background()[0] << std::endl;
  Log::test() << "Multiplicative Analysis member 1:" << anEns[0] << std::endl;

  if (params_.levels.value() != boost::none && params_.factors.value() != boost::none &&
      params_.factor.value() == boost::none) {
    Log::trace() << "MultiplicativeInflation::State inflate by level starting" << std::endl;
    this->inflateByLevel(anEns);
    Log::trace() << "MultiplicativeInflation::State inflate by level done" << std::endl;
  } else if (params_.factor.value() != boost::none && params_.levels.value() == boost::none &&
             params_.factors.value() == boost::none) {
    Log::trace() << "MultiplicativeInflation::State inflate all levels starting" << std::endl;
    double factor = *params_.factor.value();
    // calculate ensemble mean
    StateSet_ an_mean = anEns.ens_mean();

    // update analysis with MultiplicativeInflation
    for (size_t jj = 0; jj < anEns.size(); ++jj) {
      Increment_ pertTot(this->geometry(), this->vars(), anEns[jj].validTime());
      pertTot.zero();
      Increment_ pert(pertTot);

      pert.diff(anEns[jj], an_mean[0]);
      pertTot.axpy(factor, pert);

      // an_mean contains copies of anEns[0] for non-state variables
      // using zero+accumul instead of "=" ensures that only
      // state variables are modified in anEns[jj]
      anEns[jj].zero();
      anEns[jj].accumul(1.0, an_mean[0]);

      // add analysis variable perturbations
      anEns[jj] += pertTot;
  }
  Log::trace() << "MultiplicativeInflation::State inflate all levels done" << std::endl;
  } else {
    ABORT("Incorrect input settings");
  }
  Log::test() << "Multiplicative Updated Analysis member 1:" << anEns[0] << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void MultiplicativeInflation<MODEL>::inflateByLevel(IncrementSet_ & anEns) {
  IncrementSet_ an_mean = anEns.ens_mean();
  StateSet_ bgMean = this->background().ens_mean();
  Increment_ bgPert(this->geometry(), this->vars(), this->background()[0].validTime());

  std::vector<double> factors;
  this->interpolateFactors(factors, anEns[0].fieldSet().fieldSet());

  // update analysis with MultiplicativeInflation
  for (size_t jj = 0; jj < anEns.size(); ++jj) {
    bgPert.diff(this->background()[jj], bgMean[0]);
    bgPert -= an_mean[0];
    for (auto & field : anEns[jj].fieldSet().fieldSet()) {
      atlas::Field bgField = bgPert.fieldSet().fieldSet().field(field.name());
      auto view = atlas::array::make_view<double, 2>(field);
      const auto bgView = atlas::array::make_view<double, 2>(bgField);
      for (int jnode = 0; jnode < field.shape(0); ++jnode) {
        for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) *= factors[jlevel];
          double bgTerm = factors[jlevel] - 1.0;
          bgTerm *= bgView(jnode, jlevel);
          view(jnode, jlevel) += bgTerm;
        }
      }
    }
    anEns[jj].synchronizeFields();
  }
}

// -----------------------------------------------------------------------------
template<typename MODEL>
void MultiplicativeInflation<MODEL>::inflateByLevel(StateSet_ & anEns) {
  StateSet_ an_mean = anEns.ens_mean();
  Increment_ pert(this->geometry(), this->vars(), anEns[0].validTime());
  std::vector<double> factors;
  this->interpolateFactors(factors, anEns[0].fieldSet().fieldSet());

  // update analysis with MultiplicativeInflation
  for (size_t jj = 0; jj < anEns.size(); ++jj) {
    pert.diff(anEns[jj], an_mean[0]);
    for (auto & field : pert.fieldSet().fieldSet()) {
      auto view = atlas::array::make_view<double, 2>(field);
      for (int jnode = 0; jnode < field.shape(0); ++jnode) {
        for (int jlevel = 0; jlevel < field.shape(1); ++jlevel) {
          view(jnode, jlevel) *= factors[jlevel];
        }
      }
    }
    anEns[jj].zero();
    anEns[jj].accumul(1.0, an_mean[0]);

    // add analysis variable perturbations
    pert.synchronizeFields();
    anEns[jj] += pert;
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void MultiplicativeInflation<MODEL>::interpolateFactors(std::vector<double> & factors,
                                                        atlas::FieldSet & fieldset) {
  Log::trace() << "MultiplicativeInflation::interpolateFactors starting" << std::endl;
  std::vector<double> levels = *params_.levels.value();
  std::vector<double> values = *params_.factors.value();
  std::size_t nlevels = fieldset[0].shape(1);

  ASSERT(levels.size() == values.size());
  ASSERT(levels.size() <= nlevels);
  if (!std::is_sorted(levels.begin(), levels.end(), std::less_equal<double>())) {
    ABORT("Levels must be entered in ascending order with no duplication");
  }

  if (levels[0] != 0) {
    levels.insert(levels.begin(), 0);
    values.insert(values.begin(), 1.0);
  }
  if (levels.back() != nlevels - 1) {
    levels.push_back(nlevels - 1);
    values.push_back(1.0);
  }
  factors.push_back(values[0]);

  for (size_t jlevel = 1; jlevel < levels.size(); ++jlevel) {
    double iter = levels[jlevel] - levels[jlevel - 1];
    double grad = values[jlevel] - values[jlevel - 1];
    grad /= iter;
    for (int jiter = 1; jiter < iter + 1; ++jiter) {
      double factor = values[jlevel - 1] + (grad * jiter);
      factors.push_back(factor);
    }
  }
  ASSERT(factors.size() == nlevels);
  Log::trace() << "MultiplicativeInflation::factors: " << factors << std::endl;
  Log::trace() << "MultiplicativeInflation::interpolateFactors done" << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_MULTIPLICATIVEINFLATION_H_
