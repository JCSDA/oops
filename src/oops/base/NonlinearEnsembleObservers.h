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
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/Departures.h"
#include "oops/base/EnsembleObserversBase.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment4D.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/ObsEnsemble.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/StateSet.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/util/Logger.h"
#include "oops/util/printRunStats.h"
#include "oops/util/Timer.h"

namespace oops {

/// \brief Ensemble observers that compute H(x) by directly running the nonlinear observation
/// operator for each ensemble member (as opposed to \ref LinearEnsembleObservers, which
/// linearizes the observer about the ensemble mean). Used when
/// "local ensemble DA.use linear observer" is false (the default).
///
/// When constructed with the modulated-ensemble (GETKF) infrastructure active (see
/// \ref EnsembleObserversBase::modulated_), H(x) is also computed for the modulated ensemble
/// perturbations (the vertical localization eigenvectors), needed by the GETKF solvers.
template <typename MODEL, typename OBS>
class NonlinearEnsembleObservers : public EnsembleObserversBase<MODEL, OBS> {
  typedef EnsembleObserversBase<MODEL, OBS> Base_;
  typedef Departures<OBS>             Departures_;
  typedef Increment4D<MODEL>          Increment4D_;
  typedef IncrementSet<MODEL>         IncrementSet_;
  typedef ObsDataVector<OBS, int>     ObsDataInt_;
  typedef ObsEnsemble<OBS>            ObsEnsemble_;
  using typename Base_::ModelAux_;
  using typename Base_::ObsAux_;
  using typename Base_::ObsErrors_;
  using typename Base_::Observations_;
  using typename Base_::ObsSpaces_;
  using typename Base_::StateSet_;

 public:
  static const std::string classname() {return "oops::NonlinearEnsembleObservers";}

  using Base_::Base_;

  /// computes ensemble H(\p ens_xx), returns mean H(\p ens_xx), saves as hofx \p iteration
  void computeHofX(const StateSet_ & ens_xx, const StateSet_ & xxmean, size_t iteration);
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void NonlinearEnsembleObservers<MODEL, OBS>::computeHofX(const StateSet_ & ens_xx,
                                                          const StateSet_ & xxmean,
                                                          size_t iteration) {
  util::Timer timer(classname(), "computeHofX");

  ObsEnsemble_ obsens(this->obspaces_, this->nens_);
  Observations_ y_mean_xb(this->obspaces_);

  // Initialize R_ anew for each iteration
  this->R_.reset(new ObsErrors_(this->observersconf_, this->obspaces_));

  // compute and save H(x) for the mean state
  std::vector<util::DateTime> times = xxmean.times();
  util::Duration flength = times[times.size() - 1] - times[0];
  // default_tstep = 2*observation window is passed to PseudoModel as the default
  // pseudomodel time step. It is only used when StateSet has a single state, to enable
  // processing of all observations in the specified window regardless of where in
  // the time window the state is. Observations in
  // ( max(winbgn, xx.time - tstep/2); min(winend, xx.time + tstep/2) ] are
  // processed in H(x).
  util::Duration default_tstep = (this->obspaces_.windowEnd() - this->obspaces_.windowStart()) * 2;
  const ModelAux_ moderr0(this->geometry_, eckit::LocalConfiguration());
  const ObsAux_  obsaux0(this->obspaces_, this->observersconf_);

  // save QC filters, obs bias, ob errors for H(mean(Xb) (H(x_i) are saved separately)
  eckit::LocalConfiguration config;
  config.set("save hofx", false);
  config.set("save qc", true);
  config.set("save obs errors", true);
  config.set("save obs bias", true);
  config.set("iteration", std::to_string(iteration));

  this->computeHofX4D(config, xxmean, y_mean_xb, flength, default_tstep, obsaux0, moderr0,
                      *this->R_, this->qcflags_);
  y_mean_xb.save("hofx_y_mean_xb"+std::to_string(iteration));

  // QC flags and Obs errors are set to that of the H(mean(Xb))
  this->R_->save("ObsError");
  this->initializeAssimilatedMask();

  // save hofx means that hofx will be written out into ObsSpace;
  // if run computeHofX4D several times with save hofx on,
  // the hofx will be overwritten,
  // unless each time specifying iteration differently in the passed config.
  config.set("save hofx", false);
  config.set("save qc", false);
  config.set("save obs errors", false);
  config.set("save obs bias", false);

  // use temporary objects for QC flags and obs errors for ensemble members
  // to avoid overwriting the ones from the H(mean(Xb)) calculation
  std::vector<ObsDataInt_> qcflags;
  for (size_t jobs = 0; jobs < this->obspaces_.size(); ++jobs) {
    ObsDataInt_ flags(this->obspaces_[jobs], this->obspaces_[jobs].obsvariables());
    qcflags.push_back(flags);
  }
  ObsErrors_ Rmat(this->observersconf_, this->obspaces_);

  // modulation infrastructure (only used when this->modulated_ is true)
  Increment4D_ dx(this->geometry_, this->modulated_ ? this->incvars_ : ens_xx.variables(),
                 ens_xx.times());
  std::unique_ptr<IncrementSet_> Ztmp;
  if (this->modulated_) {
    std::vector<int> evMembers(this->neig_);
    std::iota(evMembers.begin(), evMembers.end(), 0);
    Ztmp = std::make_unique<IncrementSet_>(this->geometry_, this->incvars_, ens_xx.times(),
                                            xxmean.commTime(), evMembers);
  }
  Observations_ tmpObs(this->obspaces_);

  for (size_t jj = 0; jj < this->nens_; ++jj) {
    // These are recalculated for each ensemble member
    times = ens_xx.times();
    flength = times[times.size()-1] - times[0];
    default_tstep = (this->obspaces_.windowEnd() - this->obspaces_.windowStart()) * 2;

    const ModelAux_ moderr(this->geometry_, eckit::LocalConfiguration());
    const ObsAux_  obsaux(this->obspaces_, this->observersconf_);

    // Construct a single-member StateSet and populate it from ens_xx
    StateSet_ member_xx(this->geometry_, ens_xx.variables(), times, ens_xx.commTime());
    for (size_t it = 0; it < times.size(); ++it) {
      member_xx(it, 0) = ens_xx(it, jj);
    }

    this->computeHofX4D(config, member_xx, obsens[jj], flength, default_tstep,
                        obsaux, moderr, Rmat, qcflags);
    Log::test() << "H(x) for member " << jj+1 << ":" << std::endl << obsens[jj] << std::endl;
    obsens[jj].save("hofx"+std::to_string(iteration)+"_"+std::to_string(jj+1));

    if (this->modulated_) {
      // observe the modulated ensemble perturbations (vertical localization eigenvectors)
      Log::info() << " NonlinearEnsembleObservers::computeHofX starting ensemble member "
                  << jj+1 << std::endl;
      util::printRunStats("NonlinearEnsembleObservers calculate hofx");
      for (size_t it = 0; it < times.size(); ++it) {
        dx(it, 0).diff(ens_xx(it, jj), xxmean[it]);
      }
      this->vertloc_->modulateIncrement(dx, *Ztmp);
      for (size_t ieig = 0; ieig < this->neig_; ++ieig) {
        StateSet_ tmpState = xxmean;
        for (size_t it = 0; it < times.size(); ++it) {
          tmpState[it] += (*Ztmp)(it, ieig);
        }

        this->computeHofX4D(config, tmpState, tmpObs, flength, default_tstep, obsaux, moderr,
                            Rmat, qcflags);
        // mask H(x) ensemble perturbations - i.e. make sure that obs that have
        // failed QC on one ensemble member fail for all
        Departures_ tmpDep = tmpObs - y_mean_xb;
        this->updateAssimilatedMask(tmpDep);
        tmpObs.save("hofxm"+std::to_string(iteration)+"_"+std::to_string(ieig+1)+
                    "_"+std::to_string(jj+1));
      }
    }
  }

  // calculate H(x) ensemble mean
  Observations_ yb_mean(obsens.mean());
  // treat the special case of nens=1: default option: xxmean=mean(xb) then yb_mean == y_mean_xb
  // and the below is a tautology; if use control member==true, xxmean was read from the control
  // member, then using H(xxmean) (y_mean_xb) is expected by downstream applications
  if (this->nens_ == 1) {yb_mean = y_mean_xb;}
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
