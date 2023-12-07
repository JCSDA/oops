/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_COSTFCT3DVAR_H_
#define OOPS_ASSIMILATION_COSTFCT3DVAR_H_

#include <memory>

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJb3D.h"
#include "oops/assimilation/CostJo.h"
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

/// 3D-Var Cost Function
/*!
 * This class is not really necessary since it is only a special
 * case of the more general 4D-Var cost function. It is provided
 * for readability.
 */

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS> class CostFct3DVar : public CostFunction<MODEL, OBS> {
  typedef Increment<MODEL>                Increment_;
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef ControlVariable<MODEL, OBS>     CtrlVar_;
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef CostJbTotal<MODEL, OBS>         JbTotal_;
  typedef CostJo<MODEL, OBS>              CostJo_;
  typedef Geometry<MODEL>                 Geometry_;
  typedef State<MODEL>                    State_;

 public:
  CostFct3DVar(const eckit::Configuration &, const eckit::mpi::Comm &);
  CostFct3DVar(const eckit::Configuration &, const eckit::mpi::Comm &,
              std::shared_ptr<JbTotal_> &, std::shared_ptr<CostJo_> &);
  virtual ~CostFct3DVar() {}

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

  const util::TimeWindow timeWindow_;
  const eckit::mpi::Comm & comm_;
  const Geometry_ resol_;
  const Variables ctlvars_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostFct3DVar<MODEL, OBS>::CostFct3DVar(const eckit::Configuration & config,
                                       const eckit::mpi::Comm & comm)
  : CostFunction<MODEL, OBS>::CostFunction(),
    timeWindow_(config.getSubConfiguration("time window")), comm_(comm),
    resol_(eckit::LocalConfiguration(config, "geometry"), comm),
    ctlvars_(config, "analysis variables")
{
  Log::trace() << "CostFct3DVar::CostFct3DVar start" << std::endl;

  this->setupTerms(config);  // Background is read here

  Log::info() << "3DVar window: " << timeWindow_ << std::endl;
  Log::trace() << "CostFct3DVar::CostFct3DVar done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
CostFct3DVar<MODEL, OBS>::CostFct3DVar(const eckit::Configuration & config,
                                     const eckit::mpi::Comm & comm,
                                     std::shared_ptr<JbTotal_> & Jb, std::shared_ptr<CostJo_> & Jo)
  : CostFunction<MODEL, OBS>::CostFunction(),
    timeWindow_(config.getSubConfiguration("time window")), comm_(comm),
    resol_(eckit::LocalConfiguration(config, "geometry"), comm),
    ctlvars_(config, "analysis variables")
{
  Log::trace() << "CostFct3DVar::CostFct3DVar start" << std::endl;

  this->getNonConstJb() = Jb;
  this->getNonConstJo() = Jo;
  this->getJTerms().push_back(this->getNonConstJo());

  Log::info() << "3DVar window: " << timeWindow_ << std::endl;
  Log::trace() << "CostFct3DVar::CostFct3DVar done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJb3D<MODEL, OBS> * CostFct3DVar<MODEL, OBS>::newJb(const eckit::Configuration & jbConf,
                                                  const Geometry_ & resol) const {
  Log::trace() << "CostFct3DVar::newJb" << std::endl;
  return new CostJb3D<MODEL, OBS>(timeWindow_.midpoint(), jbConf, resol, ctlvars_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJo<MODEL, OBS> * CostFct3DVar<MODEL, OBS>::newJo(const eckit::Configuration & joConf) const {
  Log::trace() << "CostFct3DVar::newJo" << std::endl;
  return new CostJo<MODEL, OBS>(joConf, comm_, timeWindow_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostTermBase<MODEL, OBS> * CostFct3DVar<MODEL, OBS>::newJc(const eckit::Configuration &,
                                                           const Geometry_ &) const {
  Log::trace() << "CostFct3DVar::newJc" << std::endl;
// For now there is no Jc that can work with 3D-Var
  return nullptr;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct3DVar<MODEL, OBS>::runNL(CtrlVar_ & xx, PostProcessor<State_> & post) const {
  Log::trace() << "CostFct3DVar::runNL start" << std::endl;
  ASSERT(xx.states().is_3d());
  ASSERT(xx.state().validTime() == timeWindow_.midpoint());

  post.initialize(xx.state(), timeWindow_.midpoint(), timeWindow_.length());
  post.process(xx.state());
  post.finalize(xx.state());

  ASSERT(xx.state().validTime() == timeWindow_.midpoint());
  Log::trace() << "CostFct3DVar::runNL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFct3DVar<MODEL, OBS>::doLinearize(const Geometry_ & res, const eckit::Configuration & conf,
                                           CtrlVar_ &, CtrlVar_ &,
                                           PostProcessor<State_> & pp,
                                           PostProcessorTLAD<MODEL> & pptraj) {
  Log::trace() << "CostFct3DVar::doLinearize start" << std::endl;
  pp.enrollProcessor(new TrajectorySaver<MODEL>(conf, res, pptraj));
  Log::trace() << "CostFct3DVar::doLinearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct3DVar<MODEL, OBS>::runTLM(CtrlInc_ & dx,
                                      PostProcessorTLAD<MODEL> & cost,
                                      PostProcessor<Increment_> post,
                                      const bool) const {
  Log::trace() << "CostFct3DVar::runTLM start" << std::endl;
  ASSERT(dx.states().is_3d());
  ASSERT(dx.state().validTime() == timeWindow_.midpoint());

  cost.initializeTL(dx.state(), timeWindow_.midpoint(), timeWindow_.length());
  post.initialize(dx.state(), timeWindow_.midpoint(), timeWindow_.length());

  cost.processTL(dx.state());
  post.process(dx.state());

  cost.finalizeTL(dx.state());
  post.finalize(dx.state());

  ASSERT(dx.state().validTime() == timeWindow_.midpoint());
  Log::trace() << "CostFct3DVar::runTLM done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct3DVar<MODEL, OBS>::zeroAD(CtrlInc_ & dx) const {
  Log::trace() << "CostFct3DVar::zeroAD start" << std::endl;
  ASSERT(dx.states().is_3d());
  dx.state().zero(timeWindow_.midpoint());
  dx.modVar().zero();
  dx.obsVar().zero();
  Log::trace() << "CostFct3DVar::zeroAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct3DVar<MODEL, OBS>::runADJ(CtrlInc_ & dx,
                                      PostProcessorTLAD<MODEL> & cost,
                                      PostProcessor<Increment_> post,
                                      const bool) const {
  Log::trace() << "CostFct3DVar::runADJ start" << std::endl;
  ASSERT(dx.states().is_3d());
  ASSERT(dx.state().validTime() == timeWindow_.midpoint());

  post.initialize(dx.state(), timeWindow_.midpoint(), timeWindow_.length());
  cost.initializeAD(dx.state(), timeWindow_.midpoint(), timeWindow_.length());

  cost.processAD(dx.state());
  post.process(dx.state());

  cost.finalizeAD(dx.state());
  post.finalize(dx.state());

  ASSERT(dx.state().validTime() == timeWindow_.midpoint());

  Log::trace() << "CostFct3DVar::runADJ done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFct3DVar<MODEL, OBS>::addIncr(CtrlVar_ & xx, const CtrlInc_ & dx,
                                       PostProcessor<Increment_> &) const {
  Log::trace() << "CostFct3DVar::addIncr start" << std::endl;
  ASSERT(xx.states().is_3d());
  ASSERT(dx.states().is_3d());
  ASSERT(xx.state().validTime() == timeWindow_.midpoint());
  ASSERT(dx.state().validTime() == timeWindow_.midpoint());
  xx.state() += dx.state();
  Log::trace() << "CostFct3DVar::addIncr done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFCT3DVAR_H_
