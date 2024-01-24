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
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/State.h"
#include "oops/base/StateSet.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

namespace oops {

class RTPSParameters : public InflationParameters {
  OOPS_CONCRETE_PARAMETERS(RTPSParameters, InflationParameters);
 public:
    RequiredParameter<double> factor{"factor", this};
};

/// \brief Application for relaxation to prior spread (RTPS) inflation
///
/// \details An application updating the analysis spread using the following equation:
///
/// \f$ xa'_i <- xa'_i(((\alpha * sig_b) + sig_a(1 - \alpha))/sig_a)\f$,
///
/// where sib_b/sig_a refers to background/analysis ensemble standard deviation
/// at each grid point, \alpha is the prescribed factor, and i refers to member i.
/// Either an ensemble of states or an ensemble of increments can be inflated. If the
/// latter, each increment is added to its corresponding background before calling
/// the method for an ensemble of analysis states
///
/// See:
/// Whitaker, J. S., and T. M. Hamill, 2012: Evaluating Methods to Account for
/// System Errors in Ensemble Data Assimilation. Mon. Wea. Rev., 140, 3078–3089,
/// https://doi.org/10.1175/MWR-D-11-00276.1.

template <typename MODEL> class RTPS : public InflationBase<MODEL> {
  typedef Geometry<MODEL>                   Geometry_;
  typedef Increment<MODEL>                  Increment_;
  typedef IncrementSet<MODEL>               IncrementSet_;
  typedef State<MODEL>                      State_;
  typedef StateSet<MODEL>                   StateSet_;

 public:
  typedef RTPSParameters                    Parameters_;
  RTPS(const Parameters_ &, const Geometry_ &,
       const StateSet_ &, const Variables &);

  void doInflation(IncrementSet_ &) override;
  void doInflation(StateSet_ &) override;

 private:
  const double factor_;
};

// -----------------------------------------------------------------------------

template<typename MODEL>
RTPS<MODEL>::RTPS(const Parameters_ & params, const Geometry_ & geom,
                 const StateSet_ & bens, const Variables & vars)
  : InflationBase<MODEL>(geom, bens, vars), factor_(params.factor.value())
{
  Log::trace() << "RTPS::set up RTPS" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void RTPS<MODEL>::doInflation(IncrementSet_ & anEns) {
  StateSet_ anEnsStates(this->background());
  Log::test() << "RTPS Analysis Increment member 1:" << anEns[0] << std::endl;

  anEnsStates += anEns;

  this->doInflation(anEnsStates);

  anEns.diff(anEnsStates, this->background());

  Log::test() << "RTPS Updated Analysis Increment member 1:" << anEns[0] << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void RTPS<MODEL>::doInflation(StateSet_ & anEns) {
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

  // calculate inflation factor
  Increment_ inflation(an_stdDev[0]);
  inflation *= (1.0 - factor_);
  inflation.axpy(factor_, bg_stdDev[0]);
  inflation.fieldSet() /= an_stdDev[0].fieldSet();
  inflation.synchronizeFields();

  Log::test() << "RTPS Analysis State member 1:" << anEns[0] << std::endl;

  Increment_ pertTot(this->geometry(), this->vars(), anEns[0].validTime());
  // update analysis with RTPS
  Log::trace() << "RTPS:: update analysis with RTPS" << std::endl;
  for (size_t jj = 0; jj < anEns.size(); ++jj) {
    // calculate RTPS perturbation
    pertTot.zero();

    pertTot.diff(anEns[jj], an_mean[0]);
    pertTot.schur_product_with(inflation);

    // an_mean contains copies of anEns[0] for non-state variables
    // using zero+accumul instead of "=" ensures that only
    // state variables are modified in anEns[jj]
    anEns[jj].zero();
    anEns[jj].accumul(1.0, an_mean[0]);

    // add analysis variable perturbations
    anEns[jj] += pertTot;
  }
  Log::test() << "RTPS Updated Analysis State member 1:" << anEns[0] << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_BASE_RTPS_H_
