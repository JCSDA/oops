/*
 * (C) Crown Copyright 2023, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_COSTFCTPERT_H_
#define OOPS_ASSIMILATION_COSTFCTPERT_H_

#include <memory>

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJb3D.h"
#include "oops/assimilation/CostJo.h"
#include "oops/assimilation/CostJoPert.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/State.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

/// 3D-Var Cost Function for the Pert members of Control-Pert EDA

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS> class CostFctPert : public CostFunction<MODEL, OBS> {
  typedef Increment<MODEL>                Increment_;
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef ControlVariable<MODEL, OBS>     CtrlVar_;
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef Geometry<MODEL>                 Geometry_;
  typedef State<MODEL>                    State_;

 public:
  CostFctPert(const eckit::Configuration &, const eckit::mpi::Comm &,
              const CostFct_ & , const bool);
  virtual ~CostFctPert() {}

  void runTLM(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
              PostProcessor<Increment_>, const bool) const override;
  void runADJ(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
              PostProcessor<Increment_>, const bool) const override;
  void zeroAD(CtrlInc_ &) const override;

  void runNL(CtrlVar_ &, PostProcessor<State_>&) const override;

 private:
  void addIncr(CtrlVar_ &, const CtrlInc_ &, PostProcessor<Increment_>&) const override;

  CostJb3D<MODEL, OBS> * newJb(const eckit::Configuration &, const Geometry_ &) const override;
  CostJo<MODEL, OBS> * newJo(const eckit::Configuration &) const override;
  CostTermBase<MODEL, OBS> * newJc(const eckit::Configuration &, const Geometry_ &) const override;
  void doLinearize(const Geometry_ &, const eckit::Configuration &, CtrlVar_ &, CtrlVar_ &,
                   PostProcessor<State_> &, PostProcessorTLAD<MODEL> &) override;
  const Geometry_ & geometry() const override {return resol_;}

  util::Duration windowLength_;
  util::DateTime windowBegin_;
  util::DateTime windowEnd_;
  util::DateTime windowHalf_;
  const eckit::mpi::Comm & comm_;
  const Geometry_ resol_;
  const Variables ctlvars_;
  ControlVariable<MODEL, OBS> traj_;
  const bool shiftTime_;
  std::shared_ptr<ObserversTLAD<MODEL, OBS>> obstlad_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostFctPert<MODEL, OBS>::CostFctPert(const eckit::Configuration & config,
                                     const eckit::mpi::Comm & comm,
                                     const CostFct_ & Jlin, const bool shiftTime)
  : CostFunction<MODEL, OBS>::CostFunction(),
    windowLength_(), windowHalf_(), comm_(comm),
    resol_(eckit::LocalConfiguration(config, "geometry"), comm),
    ctlvars_(config, "analysis variables"),
    traj_(Jlin.jb().getBackground()), shiftTime_(shiftTime)
{
  Log::trace() << "CostFctPert::CostFctPert start" << std::endl;
  windowLength_ = util::Duration(config.getString("window length"));
  windowBegin_ = util::DateTime(config.getString("window begin"));
  windowEnd_ = windowBegin_ + windowLength_;
  windowHalf_ = windowBegin_ + windowLength_/2;

  obstlad_ = Jlin.jo().observersTLAD();

  this->setupTerms(config);  // The full-field ensemble member background is read here

// The difference between the ensemble member's background and the control member's background
// becomes the background for the Pert DA; this is stored as a ControlVariable object and is
// deemed to be valid at the middle of the assimilation window (as 3D-Var is used)
  ControlIncrement<MODEL, OBS> xxPertInc(Jlin.jb());
  xxPertInc.diff(this->jb().getBackground(), traj_);
  ControlVariable<MODEL, OBS> xxPert(this->jb().getBackground(), false);
  xxPert += xxPertInc;
  if (shiftTime_) {
    xxPert.state().updateTime(windowLength_/2);
    traj_.state().updateTime(windowLength_/2);
  }
  this->getNonConstJb().getBackground() = xxPert;

  Log::info() << "Pert window: begin = " << windowBegin_ << ", end = " << windowEnd_ << std::endl;
  Log::trace() << "CostFctPert::CostFctPert done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJb3D<MODEL, OBS> * CostFctPert<MODEL, OBS>::newJb(const eckit::Configuration & jbConf,
                                                      const Geometry_ & resol) const {
  Log::trace() << "CostFctPert::newJb" << std::endl;
// CostFctPert needs to use the same B matrix as the control member's cost function.
// For now, CostFctPert only supports the use of CostJb3D in the state part of the Jb term.
// This restricts the admissble options of the control member's cost function to
// CostFct3DVar, CostFctFGAT and CostFct4DVar.
  return new CostJb3D<MODEL, OBS>(shiftTime_? windowBegin_ : windowHalf_,
                                  jbConf, resol, ctlvars_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJo<MODEL, OBS> * CostFctPert<MODEL, OBS>::newJo(const eckit::Configuration & joConf) const {
  Log::trace() << "CostFctPert::newJo" << std::endl;
// CostFctPert uses a different class for the Jo term compared to the other cost functions
  return new CostJoPert<MODEL, OBS>(joConf, comm_, windowBegin_, windowEnd_, obstlad_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostTermBase<MODEL, OBS> * CostFctPert<MODEL, OBS>::newJc(const eckit::Configuration &,
                                                          const Geometry_ &) const {
  Log::trace() << "CostFctPert::newJc" << std::endl;
// For now there is no Jc that can work with 3D-Var
  return nullptr;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctPert<MODEL, OBS>::runNL(CtrlVar_ & xx, PostProcessor<State_> & post) const {
  Log::trace() << "CostFctPert::runNL start" << std::endl;
// Same as CostFct3DVar
  ASSERT(xx.states().is_3d());
  ASSERT(xx.state().validTime() == windowHalf_);

  post.initialize(xx.state(), windowHalf_, windowLength_);
  post.process(xx.state());
  post.finalize(xx.state());

  ASSERT(xx.state().validTime() == windowHalf_);
  Log::trace() << "CostFctPert::runNL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFctPert<MODEL, OBS>::doLinearize(const Geometry_ & res, const eckit::Configuration & conf,
                                          CtrlVar_ &, CtrlVar_ &,
                                          PostProcessor<State_> & pp,
                                          PostProcessorTLAD<MODEL> & pptraj) {
  Log::trace() << "CostFctPert::doLinearize start" << std::endl;
// Re-linearize the Jb term based on the linearization trajectory
  this->getNonConstJb().setPostProcTraj(traj_, conf, res, pptraj);
  pp.enrollProcessor(new TrajectorySaver<MODEL>(conf, res, pptraj));
  Log::trace() << "CostFctPert::doLinearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctPert<MODEL, OBS>::runTLM(CtrlInc_ & dx,
                                     PostProcessorTLAD<MODEL> & cost,
                                     PostProcessor<Increment_> post,
                                     const bool) const {
  Log::trace() << "CostFctPert::runTLM start" << std::endl;
// Same as CostFct3DVar
  ASSERT(dx.states().is_3d());
  ASSERT(dx.state().validTime() == windowHalf_);

  cost.initializeTL(dx.state(), windowHalf_, windowLength_);
  post.initialize(dx.state(), windowHalf_, windowLength_);

  cost.processTL(dx.state());
  post.process(dx.state());

  cost.finalizeTL(dx.state());
  post.finalize(dx.state());

  ASSERT(dx.state().validTime() == windowHalf_);
  Log::trace() << "CostFctPert::runTLM done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctPert<MODEL, OBS>::zeroAD(CtrlInc_ & dx) const {
  Log::trace() << "CostFctPert::zeroAD start" << std::endl;
// Same as CostFct3DVar
  ASSERT(dx.states().is_3d());
  dx.state().zero(windowHalf_);
  dx.modVar().zero();
  dx.obsVar().zero();
  Log::trace() << "CostFctPert::zeroAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctPert<MODEL, OBS>::runADJ(CtrlInc_ & dx,
                                     PostProcessorTLAD<MODEL> & cost,
                                     PostProcessor<Increment_> post,
                                     const bool) const {
  Log::trace() << "CostFctPert::runADJ start" << std::endl;
// Same as CostFct3DVar
  ASSERT(dx.states().is_3d());
  ASSERT(dx.state().validTime() == windowHalf_);

  post.initialize(dx.state(), windowHalf_, windowLength_);
  cost.initializeAD(dx.state(), windowHalf_, windowLength_);

  cost.processAD(dx.state());
  post.process(dx.state());

  cost.finalizeAD(dx.state());
  post.finalize(dx.state());

  ASSERT(dx.state().validTime() == windowHalf_);

  Log::trace() << "CostFctPert::runADJ done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFctPert<MODEL, OBS>::addIncr(CtrlVar_ & xx, const CtrlInc_ & dx,
                                      PostProcessor<Increment_> &) const {
  Log::trace() << "CostFctPert::addIncr start" << std::endl;
  ASSERT(xx.states().is_3d());
  ASSERT(dx.states().is_3d());
// xx and dx need not be valid at the middle of the assimilation window as long as
// their valid times agree with each other
  ASSERT(xx.state().validTime() == dx.state().validTime());
  xx.state() += dx.state();
  Log::trace() << "CostFctPert::addIncr done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFCTPERT_H_
