/*
 * (C) Copyright 2026- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <numeric>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/Departures.h"
#include "oops/base/EnsembleObserversBase.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/Increment4D.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/LinearModel.h"
#include "oops/base/ObsAuxIncrements.h"
#include "oops/base/ObsEnsemble.h"
#include "oops/base/Observations.h"
#include "oops/base/ObserversTLAD.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/StateSet.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/generic/PseudoLinearModelIncrement4D.h"
#include "oops/interface/ModelAuxIncrement.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace oops {

/// \brief Ensemble observers that compute H(x) by linearizing the observation operator about the
/// ensemble mean and applying it to ensemble perturbations (as opposed to
/// \ref NonlinearEnsembleObservers, which runs the nonlinear observer for each member). Used when
/// "local ensemble DA.use linear observer" is true.
///
/// Follows the procedure from Shlyaeva, A., & Whitaker, J. S. (2018)
/// https://doi.org/10.1029/2018MS001309 : the observation operator is linearized about the mean
/// state (cf \ref computeHofX), and applied to the ensemble perturbations from the mean
/// (cf \ref applyLinearToPerturbations). The result is added to the nonlinear observation
/// operator applied to the mean state.
///
/// When constructed with the modulated-ensemble (GETKF) infrastructure active (see
/// \ref EnsembleObserversBase::modulated_), H(x) is also (linearly) approximated for the
/// modulated ensemble perturbations (the vertical localization eigenvectors), needed by the
/// GETKF solvers.
template <typename MODEL, typename OBS>
class LinearEnsembleObservers : public EnsembleObserversBase<MODEL, OBS> {
  typedef EnsembleObserversBase<MODEL, OBS> Base_;
  typedef Departures<OBS>             Departures_;
  typedef Increment<MODEL>            Increment_;
  typedef Increment4D<MODEL>          Increment4D_;
  typedef IncrementSet<MODEL>         IncrementSet_;
  typedef LinearModel<MODEL>          LinearModel_;
  typedef ModelAuxIncrement<MODEL>    ModelAuxInc_;
  typedef ObsAuxIncrements<OBS>       ObsAuxInc_;
  typedef ObsDataVector<OBS, int>     ObsDataInt_;
  typedef ObsEnsemble<OBS>            ObsEnsemble_;
  typedef ObserversTLAD<MODEL, OBS>   ObserversTLAD_;
  typedef PseudoLinearModelIncrement4D<MODEL> PseudoLinearModel_;
  using typename Base_::ModelAux_;
  using typename Base_::ObsAux_;
  using typename Base_::ObsErrors_;
  using typename Base_::Observations_;
  using typename Base_::ObsSpaces_;
  using typename Base_::State_;
  using typename Base_::StateSet_;

 public:
  static const std::string classname() {return "oops::LinearEnsembleObservers";}

  using Base_::Base_;

  /// computes ensemble H(\p ens_xx), returns mean H(\p ens_xx), saves as hofx \p iteration
  void computeHofX(const StateSet_ & ens_xx, const StateSet_ & xxmean, size_t iteration);

 private:
  /// Runs a linear model on 4D perturbations from a state \p dx (typically an ensemble
  /// perturbation from the mean, or a modulated ensemble/vertical localization eigenvector
  /// perturbation) and applies the linearized observation operator \ref linear_hofx_ (about the
  /// trajectory set up in \ref computeHofX) to the resulting background departures, returned in
  /// \p tmpDeps.
  void applyLinearToPerturbations(const Increment4D_ & dx, const util::Duration & flength,
                                  const util::Duration & default_tstep,
                                  const ObsAuxInc_ & obsauxinc, const ModelAuxInc_ & moderrinc,
                                  const PostProcessor<Increment_> & posttl,
                                  const PostProcessorTLAD<MODEL> & posttrajtl,
                                  Departures_ & tmpDeps);

  std::unique_ptr<ObserversTLAD_> linear_hofx_;  ///< linear observer, set up in computeHofX
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void LinearEnsembleObservers<MODEL, OBS>::applyLinearToPerturbations(
    const Increment4D_ & dx, const util::Duration & flength, const util::Duration & default_tstep,
    const ObsAuxInc_ & obsauxinc, const ModelAuxInc_ & moderrinc,
    const PostProcessor<Increment_> & posttl, const PostProcessorTLAD<MODEL> & posttrajtl,
    Departures_ & tmpDeps) {
  Increment_ init_dx = dx[0];
  std::unique_ptr<PseudoLinearModel_> pseudolinearmodel =
        std::make_unique<PseudoLinearModel_>(dx, default_tstep);
  const LinearModel_ linear_model(std::move(pseudolinearmodel));
  // run linear model on the ensemble perturbation, compute linear H*dx
  linear_model.forecastTL(init_dx, moderrinc, flength, posttl, posttrajtl);
  linear_hofx_->finalizeTL(obsauxinc, tmpDeps);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void LinearEnsembleObservers<MODEL, OBS>::computeHofX(const StateSet_ & ens_xx,
                                                       const StateSet_ & xxmean,
                                                       size_t iteration) {
  util::Timer timer(classname(), "computeHofX");

  ObsEnsemble_ obsens(this->obspaces_, this->nens_);
  Observations_ y_mean_xb(this->obspaces_);

  // Initialize R_ anew for each iteration
  this->R_.reset(new ObsErrors_(this->observersconf_, this->obspaces_));

  const std::vector<util::DateTime> times = ens_xx.times();
  const util::Duration flength = times[times.size() - 1] - times[0];
  // default_tstep = 2*observation window is passed to PseudoModel as the default
  // pseudomodel time step, see NonlinearEnsembleObservers::computeHofX for more detail.
  const util::Duration default_tstep = (this->obspaces_.windowEnd()
                                        - this->obspaces_.windowStart()) * 2;
  const ModelAux_ moderr(this->geometry_, eckit::LocalConfiguration());
  const ModelAuxInc_  moderrinc(this->geometry_, eckit::LocalConfiguration());
  const ObsAux_  obsaux(this->obspaces_, this->observersconf_);
  const ObsAuxInc_  obsauxinc(this->obspaces_, this->observersconf_);

  // set up postprocessors for the linear model run on ensemble perturbations
  PostProcessor<Increment_> posttl;
  PostProcessorTLAD<MODEL> posttrajtl;

  // save QC filters, obs bias, ob errors for H(mean(Xb) (H(x_i) are saved separately)
  eckit::LocalConfiguration config;
  config.set("save hofx", false);
  config.set("save qc", true);
  config.set("save obs errors", true);
  config.set("save obs bias", true);
  config.set("iteration", std::to_string(iteration));

  // run the model forecast on the ensemble mean and linearize the observer about this trajectory.
  // save obs errors and qc flags for H(mean(Xb)). The linear observer must be set up (and its
  // trajectory-saving postprocessor enrolled) only after the nonlinear observer has been
  // initialized (i.e. inside this hook, which computeHofX4D invokes right after its own
  // nonlinear-observer initialize() and before running the forecast) to match the ordering the
  // linear observer relies on.
  PostProcessorTLAD<MODEL> posttraj;
  this->computeHofX4D(config, xxmean, y_mean_xb, flength, default_tstep, obsaux, moderr,
                      *this->R_, this->qcflags_,
                      [&](PostProcessor<State_ > & post) {
                        linear_hofx_ = std::make_unique<ObserversTLAD_>(this->obspaces_,
                                                                        this->obsconf_);
                        linear_hofx_->initializeTraj(this->geometry_, obsaux, posttraj);
                        post.enrollProcessor(new TrajectorySaver<MODEL>(
                            eckit::LocalConfiguration(), this->geometry_, posttraj));
                      });
  linear_hofx_->finalizeTraj(this->qcflags_);
  y_mean_xb.save("hofx_y_mean_xb"+std::to_string(iteration));

  // QC flags and Obs errors are set to that of the H(mean(Xb))
  this->R_->save("ObsError");
  this->initializeAssimilatedMask();

  // add linearized H(x) to the linear model postprocessor
  linear_hofx_->initializeTL(posttrajtl);

  Departures_ tmpDeps(this->obspaces_);
  Observations_ yb_mean(y_mean_xb);

  if (this->modulated_) {
    // calculate obs departures and mask out any missing departures as well as those that have
    // failed QC
    Observations_ yobs(this->obspaces_, "ObsValue");
    Departures_ omb = yobs - y_mean_xb;
    this->updateAssimilatedMask(omb);

    Increment4D_ dx(this->geometry_, this->incvars_, times);
    std::vector<int> evMembers(this->neig_);
    std::iota(evMembers.begin(), evMembers.end(), 0);
    IncrementSet_ Ztmp(this->geometry_, this->incvars_, times, xxmean.commTime(), evMembers);

    for (size_t jj = 0; jj < this->nens_; ++jj) {
      Log::info() << " LinearEnsembleObservers::computeHofX starting ensemble member "
                  << jj+1 << std::endl;
      util::printRunStats("LinearEnsembleObservers calculate hofx");
      tmpDeps.zero();
      // Setup PseudoLinearModelIncrement4D to run on ensemble perturbation
      for (size_t it = 0; it < times.size(); ++it) {
        dx(it, 0).diff(ens_xx(it, jj), xxmean[it]);
      }
      // Approximate H(x_i) around the ensemble mean using linearized model and linearized
      // observer.
      applyLinearToPerturbations(dx, flength, default_tstep, obsauxinc, moderrinc, posttl,
                                 posttrajtl, tmpDeps);
      Observations_ tmpObs(y_mean_xb);
      tmpObs += tmpDeps;
      Log::test() << "H(x) for member " << jj+1 << ":" << std::endl << tmpObs << std::endl;
      tmpObs.save("hofx"+std::to_string(iteration)+"_"+std::to_string(jj+1));

      // observe the modulated ensemble perturbations (vertical localization eigenvectors)
      this->vertloc_->modulateIncrement(dx, Ztmp);
      for (size_t ieig = 0; ieig < this->neig_; ++ieig) {
        Increment4D_ z_eig(this->geometry_, this->incvars_, times);
        for (size_t it = 0; it < times.size(); ++it) {
          z_eig[it] = Ztmp(it, ieig);
        }
        applyLinearToPerturbations(z_eig, flength, default_tstep, obsauxinc, moderrinc,
                                   posttl, posttrajtl, tmpDeps);
        Observations_ tmpObsEig(y_mean_xb);
        tmpObsEig += tmpDeps;
        tmpObsEig.save("hofxm"+std::to_string(iteration)+"_"+std::to_string(ieig+1)+
                       "_"+std::to_string(jj+1));
      }
    }
    // for linear H, the ensemble mean H(x) is equal to H(mean(Xb)), computed above
  } else {
    for (size_t jj = 0; jj < this->nens_; ++jj) {
      // Setup PseudoLinearModelIncrement4D to run on ensemble perturbation
      Increment4D_ dx(this->geometry_, ens_xx.variables(), times);
      for (size_t it = 0; it < times.size(); ++it) {
        dx(it, 0).diff(ens_xx(it, jj), xxmean[it]);
      }
      // Approximate H(x_i) (obsens[jj]) around the ensemble mean using linearized model and
      // linearized observer. Firstly, apply the linearized obs operator to this ensemble
      // member's background perturbation from the ensemble mean.
      applyLinearToPerturbations(dx, flength, default_tstep, obsauxinc, moderrinc, posttl,
                                 posttrajtl, tmpDeps);
      // Secondly, add this to the ensemble mean in observation space (calculated with the
      // nonlinear obs operator) giving the approximate H(x_i)
      obsens[jj] = y_mean_xb;
      obsens[jj] += tmpDeps;
      Log::test() << "H(x) for member " << jj+1 << ":" << std::endl << obsens[jj] << std::endl;
      obsens[jj].save("hofx"+std::to_string(iteration)+"_"+std::to_string(jj+1));
    }
    // calculate H(x) ensemble mean
    yb_mean = obsens.mean();
    // treat the special case of nens=1: default option: xxmean=mean(xb) then yb_mean == y_mean_xb
    // and the below is a tautology; if use control member==true, xxmean was read from the control
    // member, then using H(xxmean) (y_mean_xb) is expected by downstream applications
    if (this->nens_ == 1) {yb_mean = y_mean_xb;}
  }

  this->maskWithPriorQcFlags(yb_mean);
  const std::string type = iteration == 0 ? "background" : "analysis";
  const std::string group = iteration == 0 ? "ombg" : "oman";
  Log::test() << "H(x) ensemble " << type << " mean: " << std::endl << yb_mean << std::endl;
  Observations_ yobs(this->obspaces_, "ObsValue");
  Departures_ diff(yobs - yb_mean);
  this->maskWithPriorQcFlags(diff);

  diff.save(group);
  Log::test() << type << " y - H(x): " << std::endl << diff << std::endl;
  // display overall background/analysis RMS stats
  Log::test() << group << " RMS: " << diff.rms() << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops
