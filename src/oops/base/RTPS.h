/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_RTPS_H_
#define OOPS_BASE_RTPS_H_

#include <memory>
#include <string>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/InflationBase.h"
#include "oops/base/State.h"
#include "oops/base/StateSet.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

namespace oops {

/// \brief Application for relaxation to prior spread (RTPS) inflation
///
/// \details An application updating the analysis spread.
///
/// We define \gamma such that:
///
/// \f$ \gamma = 1 - \alpha + \alpha(sig_b/sig_a) \f$,
///
/// where sib_b/sig_a refers to background/analysis ensemble standard deviation
/// at each grid point and \alpha is the prescribed inflation factor.
///
/// For an ensemble of analysis states, the update is done using the following equation:
///
/// \f$ xa'_i <- \gamma * xa'_i \f$,
///
/// where xa'_i is the analysis perturbation from the mean.
///
/// For an ensemble of increments:
///
/// \f$ dxa_i <- (\gamma * dxa_i) + (\gamma - 1) * (xb'_i - dxa_mean)\f$,
///
/// where where dxa_i refers to the analysis increment for member i,
/// dxa_mean is the mean field of the ensemble of increments and
/// xb'_i is the background perturbation from the mean.
///
/// See section 3:
/// Inverarity, G.W., Tennant, W.J., Anton, L., Bowler, N.E., Clayton, A.M., Jardak, M.,
/// et al. (2023) Met Office MOGREPS-G initialisation using an ensemble of hybrid
/// four-dimensional ensemble variational (En-4DEnVar) data assimilations.
/// Quarterly Journal of the Royal Meteorological Society, 149(753), 1138–1164.
/// Available from: https://doi.org/10.1002/qj.4431

template <typename MODEL> class RTPS : public InflationBase<MODEL> {
  typedef Geometry<MODEL>                   Geometry_;
  typedef Increment<MODEL>                  Increment_;
  typedef IncrementSet<MODEL>               IncrementSet_;
  typedef State<MODEL>                      State_;
  typedef StateSet<MODEL>                   StateSet_;

 public:
  RTPS(const eckit::Configuration &, const Geometry_ &,
       const StateSet_ &, const Variables &);

  void doInflation(IncrementSet_ &) override;
  void doInflation(StateSet_ &) override;
  Increment_ computeMultiplier(const Increment_ &, const Increment_ &);

 private:
  const double factor_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
RTPS<MODEL>::RTPS(const eckit::Configuration & conf, const Geometry_ & geom,
                  const StateSet_ & bens, const Variables & vars)
  : InflationBase<MODEL>(geom, bens, vars), factor_(conf.getDouble("factor"))
{
  Log::trace() << "RTPS::set up RTPS" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void RTPS<MODEL>::doInflation(IncrementSet_ & anEns) {
  Log::trace() << "RTPS::doInflation (Increment) start" << std::endl;
  Log::test() << "RTPS Analysis Increment member 1:" << anEns[0] << std::endl;

  IncrementSet_ anEnsOriginal(anEns);

  StateSet_ bgStates(this->background());
  StateSet_ anEnsStates(this->background());
  anEnsStates += anEns;

  // calculate ensemble means
  IncrementSet_ an_mean = anEns.ens_mean();
  StateSet_ bg_mean = bgStates.ens_mean();

  // calculate ensemble standard deviations
  IncrementSet_ an_inc(this->geometry(), this->vars(), anEnsStates);
  IncrementSet_ an_stdDev(this->geometry(), this->vars(), anEns.times(), anEns.commTime());
  an_stdDev = an_inc.ens_stddev();
  IncrementSet_ bg_inc(this->geometry(), this->vars(), bgStates);
  IncrementSet_ bg_stdDev(this->geometry(), this->vars(), this->background().times(),
                          this->background().commTime());
  bg_stdDev = bg_inc.ens_stddev();

  Increment_ multiplier(anEns[0]);
  multiplier = this->computeMultiplier(an_stdDev[0], bg_stdDev[0]);
  Increment_ multiplierMinusOne(multiplier);
  multiplierMinusOne.ones();
  multiplierMinusOne *= -1;
  multiplierMinusOne += multiplier;

  Increment_ secondTerm(this->geometry(), this->vars(), anEns[0].validTime());

  for (size_t jj = 0; jj < anEns.size(); ++jj) {
    // calculate first term (gamma * analysis increment)
    anEns[jj].schur_product_with(multiplier);

    // calculate second term ((gamma - 1)*(bg perturbation minus analysis mean))
    secondTerm.diff(bgStates[jj], bg_mean[0]);
    secondTerm -= an_mean[0];
    secondTerm.schur_product_with(multiplierMinusOne);

    // add the terms together
    anEns[jj] += secondTerm;
  }

  Log::test() << "RTPS Updated Analysis Increment member 1:" << anEns[0] << std::endl;
  Log::trace() << "RTPS::doInflation (Increment) done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void RTPS<MODEL>::doInflation(StateSet_ & anEns) {
  Log::trace() << "RTPS::doInflation (State) start" << std::endl;
  // calculate ensemble mean
  StateSet_ an_mean = anEns.ens_mean();

  // calculate ensemble standard deviations
  IncrementSet_ an_inc(this->geometry(), this->vars(), anEns);
  IncrementSet_ an_stdDev(this->geometry(), this->vars(), anEns.times(), anEns.commTime());
  an_stdDev = an_inc.ens_stddev();
  IncrementSet_ bg_inc(this->geometry(), this->vars(), this->background());
  IncrementSet_ bg_stdDev(this->geometry(), this->vars(), this->background().times(),
                          this->background().commTime());
  bg_stdDev = bg_inc.ens_stddev();

  Log::test() << "RTPS Analysis State member 1:" << anEns[0] << std::endl;

  Increment_ pertTot(this->geometry(), this->vars(), anEns[0].validTime());
  Increment_ multiplier(pertTot);
  multiplier = this->computeMultiplier(an_stdDev[0], bg_stdDev[0]);
  // update analysis with RTPS
  Log::trace() << "RTPS:: update analysis with RTPS" << std::endl;
  for (size_t jj = 0; jj < anEns.size(); ++jj) {
    // calculate RTPS perturbation
    pertTot.zero();

    pertTot.diff(anEns[jj], an_mean[0]);
    pertTot.schur_product_with(multiplier);

    // an_mean contains copies of anEns[0] for non-state variables
    // using zero+accumul instead of "=" ensures that only
    // state variables are modified in anEns[jj]
    anEns[jj].zero();
    anEns[jj].accumul(1.0, an_mean[0]);

    // add analysis variable perturbations
    anEns[jj] += pertTot;
  }
  Log::test() << "RTPS Updated Analysis State member 1:" << anEns[0] << std::endl;
  Log::trace() << "RTPS::doInflation (State) done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> RTPS<MODEL>::computeMultiplier(const Increment_ & an_stdDev,
                                                        const Increment_ & bg_stdDev) {
  Log::trace() << "RTPS::computeMultiplier start" << std::endl;
  Increment_ inflation(an_stdDev);
  inflation *= (1.0 - factor_);
  inflation.axpy(factor_, bg_stdDev);
  inflation.fieldSet() /= an_stdDev.fieldSet();
  inflation.synchronizeFields();

  Log::trace() << "RTPS::computeMultiplier done" << std::endl;
  return inflation;
}

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_BASE_RTPS_H_
