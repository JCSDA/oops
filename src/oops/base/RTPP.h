/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_RTPP_H_
#define OOPS_BASE_RTPP_H_

#include <memory>
#include <string>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/InflationBase.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/State.h"
#include "oops/base/StateSet.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

namespace oops {

class RTPPParameters : public InflationParameters {
  OOPS_CONCRETE_PARAMETERS(RTPPParameters, InflationParameters);

 public:
    RequiredParameter<double> factor{"factor", this};
};

/// \brief Application for relaxation to prior perturbation (RTPP) inflation
///
/// \details An application updating the analysis perturbation.
/// For an ensemble of states it uses the following equation:
///
/// \f$ xa'_i <- (\alpha * xb'_i) + (1 - \alpha) * xa'_i\f$,
///
/// where xb'_i/xa'_i refers to background/analysis perturbation from the mean
/// for member i, and \alpha is the prescribed factor.
///
/// For an ensemble of analysis increments:
///
/// \f$ dxa_i <- (\alpha * dxa_mean) + (1 - \alpha) * dxa_i\f$,
///
/// where dxa_i refers to the analysis increment for member i,
/// dxa_mean is the mean field of the ensemble of increments,
/// and \alpha is the prescribed factor.
///
/// See section 3:
/// Inverarity, G.W., Tennant, W.J., Anton, L., Bowler, N.E., Clayton, A.M., Jardak, M.,
/// et al. (2023) Met Office MOGREPS-G initialisation using an ensemble of hybrid
/// four-dimensional ensemble variational (En-4DEnVar) data assimilations.
/// Quarterly Journal of the Royal Meteorological Society, 149(753), 1138–1164.
/// Available from: https://doi.org/10.1002/qj.4431

template <typename MODEL> class RTPP : public InflationBase<MODEL> {
  typedef Geometry<MODEL>                   Geometry_;
  typedef Increment<MODEL>                  Increment_;
  typedef IncrementSet<MODEL>               IncrementSet_;
  typedef State<MODEL>                      State_;
  typedef StateSet<MODEL>                   StateSet_;

 public:
  typedef RTPPParameters                    Parameters_;
  RTPP(const Parameters_ &, const Geometry_ &,
       const StateSet_ &, const Variables &);

  void doInflation(IncrementSet_ &) override;
  void doInflation(StateSet_ &) override;

 private:
  const double factor_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
RTPP<MODEL>::RTPP(const Parameters_ & params, const Geometry_ & geom,
                 const StateSet_ & bens, const Variables & vars)
  : InflationBase<MODEL>(geom, bens, vars), factor_(params.factor.value())
{
  Log::trace() << "RTPP::set up RTPP" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void RTPP<MODEL>::doInflation(IncrementSet_ & anEns) {
  IncrementSet_ an_mean = anEns.ens_mean();
  Log::test() << "RTPP Analysis Increment member 1:" << anEns[0] << std::endl;

  for (size_t jj = 0; jj < anEns.size(); ++jj) {
    // calculate RTPP perturbations
    anEns[jj] *= (1.0 - factor_);
    anEns[jj].axpy(factor_, an_mean[0]);
  }
  Log::test() << "RTPP Updated Analysis Increment member 1:" << anEns[0] << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void RTPP<MODEL>::doInflation(StateSet_ & anEns) {
  // calculate ensemble means
  StateSet_ bg_mean = this->background().ens_mean();
  StateSet_ an_mean = anEns.ens_mean();

  Log::test() << "RTPP Background member 1:" << this->background()[0] << std::endl;
  Log::test() << "RTPP Analysis member 1:" << anEns[0] << std::endl;

  // update analysis with RTPP
  for (size_t jj = 0; jj < anEns.size(); ++jj) {
    // calculate RTPP perturbation
    Increment_ pertTot(this->geometry(), this->vars(), anEns[jj].validTime());
    pertTot.zero();

    Increment_ pert(pertTot);

    pert.diff(this->background()[jj], bg_mean[0]);
    pertTot.axpy(factor_, pert);

    pert.diff(anEns[jj], an_mean[0]);
    pertTot.axpy((1.0 - factor_), pert);

    // an_mean contains copies of anEns[0] for non-state variables
    // using zero+accumul instead of "=" ensures that only
    // state variables are modified in anEns[jj]
    anEns[jj].zero();
    anEns[jj].accumul(1.0, an_mean[0]);

    // add analysis variable perturbations
    anEns[jj] += pertTot;
  }
  Log::test() << "RTPP Updated Analysis member 1:" << anEns[0] << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_BASE_RTPP_H_
