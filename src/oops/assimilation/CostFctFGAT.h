/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJb3D.h"
#include "oops/assimilation/CostJo.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/Model.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/State.h"
#include "oops/base/StateSaver.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/base/Variables.h"
#include "oops/interface/VariableChange.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

/// 3D-FGAT Cost Function

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS> class CostFctFGAT : public CostFunction<MODEL, OBS> {
  typedef Increment<MODEL>                Increment_;
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef ControlVariable<MODEL, OBS>     CtrlVar_;
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef Geometry<MODEL>                 Geometry_;
  typedef State<MODEL>                    State_;
  typedef Model<MODEL>                    Model_;
  typedef VariableChange<MODEL>           VarCha_;

 public:
  CostFctFGAT(const eckit::Configuration &, const eckit::mpi::Comm &);
  virtual ~CostFctFGAT() {}

  void runTLM(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
              PostProcessor<Increment_>, const bool) const override;
  void runADJ(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
              PostProcessor<Increment_>, const bool) const override;
  void zeroAD(CtrlInc_ &) const override;

  void runNL(CtrlVar_ &, PostProcessor<State_>&) const override;

 protected:
  const Geometry_ & geometry() const override {return resol_;}

 private:
  void addIncr(CtrlVar_ &, const CtrlInc_ &, PostProcessor<Increment_>&) const override;

  CostJb3D<MODEL, OBS> * newJb(const eckit::Configuration &, const Geometry_ &) const override;
  CostJo<MODEL, OBS>       * newJo(const eckit::Configuration &) const override;
  CostTermBase<MODEL, OBS> * newJc(const eckit::Configuration &, const Geometry_ &) const override;
  void doLinearize(const Geometry_ &, const eckit::Configuration &, CtrlVar_ &, CtrlVar_ &,
                   PostProcessor<State_> &, PostProcessorTLAD<MODEL> &) override;
  void finishLinearize() override;

  util::Duration windowLength_;
  util::DateTime windowBegin_;
  util::DateTime windowEnd_;
  util::DateTime windowHalf_;
  const eckit::mpi::Comm & comm_;
  const Geometry_ resol_;
  Model_ model_;
  const Variables ctlvars_;
  std::unique_ptr<VarCha_> an2model_;
  mutable bool fgat_;
  State_ * hackBG_;
  State_ * hackFG_;
  std::shared_ptr<StateSaver<State_>> saver_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostFctFGAT<MODEL, OBS>::CostFctFGAT(const eckit::Configuration & config,
                                     const eckit::mpi::Comm & comm)
  : CostFunction<MODEL, OBS>::CostFunction(), comm_(comm),
    resol_(eckit::LocalConfiguration(config, "geometry"), comm),
    model_(resol_, eckit::LocalConfiguration(config, "model")),
    ctlvars_(config, "analysis variables"),
    fgat_(false), hackBG_(nullptr), hackFG_(nullptr), saver_()
{
  Log::trace() << "CostFctFGAT::CostFctFGAT start" << std::endl;
  windowLength_ = util::Duration(config.getString("window length"));
  windowBegin_ = util::DateTime(config.getString("window begin"));
  windowEnd_ = windowBegin_ + windowLength_;
  windowHalf_ = windowBegin_ + windowLength_/2;
  an2model_ = std::make_unique<VarCha_>(config.getSubConfiguration("variable change"), resol_);

  this->setupTerms(config);  // Background is read here

  Log::info() << "FGAT window: begin = " << windowBegin_ << ", end = " << windowEnd_ << std::endl;
  Log::trace() << "CostFctFGAT::CostFctFGAT done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJb3D<MODEL, OBS> * CostFctFGAT<MODEL, OBS>::newJb(const eckit::Configuration & jbConf,
                                                 const Geometry_ & resol) const {
  Log::trace() << "CostFctFGAT::newJb" << std::endl;
  CostJb3D<MODEL, OBS> * jb = new CostJb3D<MODEL, OBS>(windowBegin_, jbConf, resol, ctlvars_);
  return jb;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJo<MODEL, OBS> * CostFctFGAT<MODEL, OBS>::newJo(const eckit::Configuration & joConf) const {
  Log::trace() << "CostFctFGAT::newJo" << std::endl;
  return new CostJo<MODEL, OBS>(joConf, comm_, windowBegin_, windowEnd_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostTermBase<MODEL, OBS> * CostFctFGAT<MODEL, OBS>::newJc(const eckit::Configuration & jcConf,
                                                          const Geometry_ &) const {
  Log::trace() << "CostFctFGAT::newJc" << std::endl;
// For now there is no Jc that can work with 3D-Var
  return nullptr;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctFGAT<MODEL, OBS>::runNL(CtrlVar_ & xx, PostProcessor<State_> & post) const {
  Log::trace() << "CostFctFGAT::runNL start" << std::endl;
  ASSERT(xx.states().is_3d());

  if (fgat_) {
    ASSERT(xx.state().validTime() == windowBegin_);

    model_.forecast(xx.state(), xx.modVar(), windowLength_, post);

    ASSERT(xx.state().validTime() == windowEnd_);
  } else {
    ASSERT(xx.state().validTime() == windowHalf_);

    post.initialize(xx.state(), windowHalf_, windowLength_);
    post.process(xx.state());
    post.finalize(xx.state());

    ASSERT(xx.state().validTime() == windowHalf_);
  }

  Log::trace() << "CostFctFGAT::runNL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFctFGAT<MODEL, OBS>::doLinearize(const Geometry_ & res, const eckit::Configuration & conf,
                                          CtrlVar_ & bg, CtrlVar_ & fg,
                                          PostProcessor<State_> & pp,
                                          PostProcessorTLAD<MODEL> & pptraj) {
  Log::trace() << "CostFctFGAT::doLinearize start" << std::endl;
  ASSERT(bg.states().is_3d());
  ASSERT(fg.states().is_3d());
  fgat_ = (conf.getInt("iteration") == 0);
  pp.enrollProcessor(new TrajectorySaver<MODEL>(conf, res, pptraj));
  hackBG_ = &bg.state();
  hackFG_ = &fg.state();

  std::vector<std::string> antime = {windowHalf_.toString()};
  eckit::LocalConfiguration halfwin;
  halfwin.set("times", antime);
  saver_.reset(new StateSaver<State_>(halfwin));
  pp.enrollProcessor(saver_);

  Log::trace() << "CostFctFGAT::doLinearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFctFGAT<MODEL, OBS>::finishLinearize() {
  Log::trace() << "CostFctFGAT::finishLinearize start" << std::endl;
  ASSERT(saver_->getState().validTime() == windowHalf_);
  Variables anvars(hackBG_->variables());
  an2model_->changeVarInverse(saver_->getState(), anvars);
  *hackBG_ = saver_->getState();
  *hackFG_ = saver_->getState();
  fgat_ = false;
  Log::trace() << "CostFctFGAT::finishLinearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctFGAT<MODEL, OBS>::runTLM(CtrlInc_ & dx,
                                     PostProcessorTLAD<MODEL> & cost,
                                     PostProcessor<Increment_> post,
                                     const bool) const {
  Log::trace() << "CostFctFGAT::runTLM start" << std::endl;
  ASSERT(dx.states().is_3d());
  ASSERT(dx.state().validTime() == windowHalf_);

  cost.initializeTL(dx.state(), windowHalf_, windowLength_);
  post.initialize(dx.state(), windowHalf_, windowLength_);

  cost.processTL(dx.state());
  post.process(dx.state());

  cost.finalizeTL(dx.state());
  post.finalize(dx.state());

  ASSERT(dx.state().validTime() == windowHalf_);
  Log::trace() << "CostFctFGAT::runTLM done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctFGAT<MODEL, OBS>::zeroAD(CtrlInc_ & dx) const {
  Log::trace() << "CostFctFGAT::zeroAD start" << std::endl;
  ASSERT(dx.states().is_3d());
  dx.state().zero(windowHalf_);
  dx.modVar().zero();
  dx.obsVar().zero();
  Log::trace() << "CostFctFGAT::zeroAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctFGAT<MODEL, OBS>::runADJ(CtrlInc_ & dx,
                                     PostProcessorTLAD<MODEL> & cost,
                                     PostProcessor<Increment_> post,
                                     const bool) const {
  Log::trace() << "CostFctFGAT::runADJ start" << std::endl;
  ASSERT(dx.states().is_3d());
  ASSERT(dx.state().validTime() == windowHalf_);

  post.initialize(dx.state(), windowHalf_, windowLength_);
  cost.initializeAD(dx.state(), windowHalf_, windowLength_);

  cost.processAD(dx.state());
  post.process(dx.state());

  cost.finalizeAD(dx.state());
  post.finalize(dx.state());

  ASSERT(dx.state().validTime() == windowHalf_);

  Log::trace() << "CostFctFGAT::runADJ done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFctFGAT<MODEL, OBS>::addIncr(CtrlVar_ & xx, const CtrlInc_ & dx,
                                      PostProcessor<Increment_> &) const {
  Log::trace() << "CostFctFGAT::addIncr start" << std::endl;
  ASSERT(xx.states().is_3d());
  ASSERT(dx.states().is_3d());
  ASSERT(xx.state().validTime() == windowHalf_);
  ASSERT(dx.state().validTime() == windowHalf_);
  xx.state() += dx.state();
  Log::trace() << "CostFctFGAT::addIncr done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

