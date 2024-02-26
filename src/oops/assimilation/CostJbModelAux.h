/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>

#include "eckit/config/LocalConfiguration.h"

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/base/Geometry.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/State.h"
#include "oops/interface/ModelAuxCovariance.h"

namespace oops {

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS> class CostJbModelAux {
  typedef ControlIncrement<MODEL, OBS>  CtrlInc_;
  typedef ControlVariable<MODEL, OBS>   CtrlVar_;
  typedef Geometry<MODEL>               Geometry_;
  typedef ModelAuxControl<MODEL>        ModelAuxControl_;
  typedef ModelAuxCovariance<MODEL>     ModelAuxCovariance_;
  typedef ModelAuxIncrement<MODEL>      ModelAuxIncrement_;
  typedef State<MODEL>                  State_;
  typedef PostProcessor<State_>         PostProc_;
  typedef PostProcessorTLAD<MODEL>      PostProcTLAD_;

 public:
  CostJbModelAux(const eckit::Configuration &, const Geometry_ &);
  ~CostJbModelAux() {}

  void setPostProcTraj(const CtrlVar_ &, const eckit::Configuration &,
                       const Geometry_ &, PostProcTLAD_ &);
  void computeCostTraj();

  void multiplyCovar(const CtrlInc_ &, CtrlInc_ &) const;
  void multiplyCoInv(const CtrlInc_ &, CtrlInc_ &) const;

  void randomize(CtrlInc_ &) const;

  std::shared_ptr<ModelAuxControl_> background() const {return bg_;}
  const Geometry_ & geometry() const {ASSERT(innerRes_); return *innerRes_;}
  const eckit::Configuration & config() const {return conf_;}

 private:
  ModelAuxCovariance_ B_;
  std::shared_ptr<ModelAuxControl_> bg_;
  const Geometry_ * innerRes_;
  const ModelAuxControl_ * traj_;
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
CostJbModelAux<MODEL, OBS>::CostJbModelAux(const eckit::Configuration & conf,
                                           const Geometry_ & resol)
  : B_(conf.getSubConfiguration("model aux error"), resol),
    bg_(new ModelAuxControl_(resol, conf.getSubConfiguration("model aux control"))),
    conf_(conf.getSubConfiguration("model aux error"))
{
  Log::trace() << "CostJbModelAux contructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbModelAux<MODEL, OBS>::setPostProcTraj(const CtrlVar_ & traj,
                                                 const eckit::Configuration &,
                                                 const Geometry_ & resol, PostProcTLAD_ &) {
  traj_ = &traj.modVar();
  innerRes_ = &resol;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbModelAux<MODEL, OBS>::computeCostTraj() {
  ASSERT(traj_);
  B_.linearize(*traj_, *innerRes_);
  traj_ = nullptr;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbModelAux<MODEL, OBS>::multiplyCovar(const CtrlInc_ & dxi, CtrlInc_ & dxo) const {
  B_.multiply(dxi.modVar(), dxo.modVar());
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbModelAux<MODEL, OBS>::multiplyCoInv(const CtrlInc_ & dxi, CtrlInc_ & dxo) const {
  B_.inverseMultiply(dxi.modVar(), dxo.modVar());
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbModelAux<MODEL, OBS>::randomize(CtrlInc_ & dx) const {
  B_.randomize(dx.modVar());
}

// -----------------------------------------------------------------------------

}  // namespace oops
