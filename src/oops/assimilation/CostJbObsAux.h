/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/base/Geometry.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsAuxCovariances.h"
#include "oops/base/ObsAuxIncrements.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/State.h"

namespace oops {

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS> class CostJbObsAux {
  typedef ControlIncrement<MODEL, OBS>  CtrlInc_;
  typedef ControlVariable<MODEL, OBS>   CtrlVar_;
  typedef Geometry<MODEL>               Geometry_;
  typedef ObsAuxControls<OBS>           ObsAuxControls_;
  typedef ObsAuxCovariances<OBS>        ObsAuxCovariances_;
  typedef ObsAuxIncrements<OBS>         ObsAuxIncrements_;
  typedef ObsSpaces<OBS>                ObsSpaces_;
  typedef State<MODEL>                  State_;
  typedef PostProcessor<State_>         PostProc_;
  typedef PostProcessorTLAD<MODEL>      PostProcTLAD_;

 public:
  CostJbObsAux(const ObsSpaces_ & odb, const eckit::Configuration &);
  ~CostJbObsAux() {}

  void setPostProcTraj(const CtrlVar_ &, const eckit::Configuration &,
                       const Geometry_ &, PostProcTLAD_ &);
  void computeCostTraj();

  void multiplyCovar(const CtrlInc_ &, CtrlInc_ &) const;
  void multiplyCoInv(const CtrlInc_ &, CtrlInc_ &) const;

  void randomize(CtrlInc_ &) const;

  std::shared_ptr<ObsAuxControls_> background() const {return bg_;}
  const ObsAuxCovariances_ & covariance() const {return B_;}
  const ObsSpaces_ & obspaces() const {return B_.obspaces();}
  const eckit::Configuration & config() const {return B_.config();}

 private:
  ObsAuxCovariances_ B_;
  std::shared_ptr<ObsAuxControls_> bg_;
  eckit::LocalConfiguration innerConf_;
  const ObsAuxControls_ * traj_;
  std::vector<bool> coldStart_;
  bool coldStartSynced_;
};

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
CostJbObsAux<MODEL, OBS>::CostJbObsAux(const ObsSpaces_ & odb, const eckit::Configuration & conf)
  : B_(odb, conf), bg_(new ObsAuxControls_(odb, conf)), coldStartSynced_(false)
{
  // Record which obs spaces asked for a VarBC cold start, so that the background can be brought
  // into line with the cold-started first guess once it has been computed (see computeCostTraj).
  const std::vector<eckit::LocalConfiguration> obsconf = conf.getSubConfigurations();
  coldStart_.reserve(obsconf.size());
  for (std::size_t jobs = 0; jobs < obsconf.size(); ++jobs) {
    const eckit::LocalConfiguration biasconf = obsconf[jobs].getSubConfiguration("obs bias");
    coldStart_.push_back(biasconf.getBool("cold start.enable", false));
  }
  Log::trace() << "CostJbObsAux contructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbObsAux<MODEL, OBS>::setPostProcTraj(const CtrlVar_ & traj,
                                               const eckit::Configuration & conf,
                                               const Geometry_ & resol, PostProcTLAD_ &) {
  traj_ = &traj.obsVar();
  innerConf_ = eckit::LocalConfiguration(conf);
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbObsAux<MODEL, OBS>::computeCostTraj() {
  ASSERT(traj_);
  // A VarBC cold start replaces zero coefficients with an estimate derived from the uncorrected
  // departures, and that estimate is what the analysis should be drawn towards. Adopt it as the
  // background before Jb is evaluated, otherwise Jb would penalize the cold start itself and
  // pull the coefficients straight back to zero. Done once, on the first outer loop, since the
  // cold start only ever runs during the first H(x).
  if (!coldStartSynced_) {
    for (std::size_t jj = 0; jj < bg_->size() && jj < coldStart_.size(); ++jj) {
      if (coldStart_[jj]) {
        (*bg_)[jj] = (*traj_)[jj];
        Log::info() << "CostJbObsAux: adopted cold-started VarBC coefficients as background "
                    << "for obs space " << jj << std::endl;
      }
    }
    coldStartSynced_ = true;
  }
  B_.linearize(*traj_, innerConf_);
  traj_ = nullptr;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbObsAux<MODEL, OBS>::multiplyCovar(const CtrlInc_ & dxi, CtrlInc_ & dxo) const {
  B_.multiply(dxi.obsVar(), dxo.obsVar());
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbObsAux<MODEL, OBS>::multiplyCoInv(const CtrlInc_ & dxi, CtrlInc_ & dxo) const {
  B_.inverseMultiply(dxi.obsVar(), dxo.obsVar());
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbObsAux<MODEL, OBS>::randomize(CtrlInc_ & dx) const {
  B_.randomize(dx.obsVar());
}

// -----------------------------------------------------------------------------

}  // namespace oops
