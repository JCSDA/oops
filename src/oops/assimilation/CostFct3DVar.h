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
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJb3D.h"
#include "oops/assimilation/CostJo.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/State.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/base/Variables.h"
#include "oops/generic/ModelBase.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

/// Parameters for the 3D-Var cost function
template <typename MODEL, typename OBS>
class CostFct3DVarParameters : public CostFunctionParametersBase<MODEL, OBS> {
  // This typedef prevents the macro below from choking on the 2 args of the templated type
  typedef CostFunctionParametersBase<MODEL, OBS> CostFuntionParametersBase_;
  OOPS_CONCRETE_PARAMETERS(CostFct3DVarParameters, CostFuntionParametersBase_);

 public:
  typedef typename State<MODEL>::Parameters_           StateParameters_;
  typedef ModelSpaceCovarianceParametersWrapper<MODEL> CovarianceParameters_;

  // options for Jb term
  RequiredParameter<StateParameters_> background{"background", "background state", this};
  RequiredParameter<CovarianceParameters_> backgroundError{"background error", "background error",
      this};
};

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
  typedef Geometry<MODEL>                 Geometry_;
  typedef State<MODEL>                    State_;

 public:
  typedef CostFct3DVarParameters<MODEL, OBS> Parameters_;

  CostFct3DVar(const Parameters_ &, const eckit::mpi::Comm &);
  virtual ~CostFct3DVar() {}

  void runTLM(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
              PostProcessor<Increment_>, const bool) const override;
  void runADJ(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
              PostProcessor<Increment_>, const bool) const override;
  void zeroAD(CtrlInc_ &) const override;

  void runNL(CtrlVar_ &, PostProcessor<State_>&) const override;

 private:
  void addIncr(CtrlVar_ &, const CtrlInc_ &, PostProcessor<Increment_>&) const override;

  CostJb3D<MODEL>     * newJb(const eckit::Configuration &, const Geometry_ &,
                              const CtrlVar_ &) const override;
  CostJo<MODEL, OBS>       * newJo(const ObserversParameters<MODEL, OBS> &) const override;
  CostTermBase<MODEL, OBS> * newJc(const eckit::Configuration &, const Geometry_ &) const override;
  void doLinearize(const Geometry_ &, const eckit::Configuration &,
                   const CtrlVar_ &, const CtrlVar_ &,
                   PostProcessor<State_> &, PostProcessorTLAD<MODEL> &) override;
  const Geometry_ & geometry() const override {return resol_;}

  util::Duration windowLength_;
  util::DateTime windowBegin_;
  util::DateTime windowEnd_;
  util::DateTime windowHalf_;
  const eckit::mpi::Comm & comm_;
  const Geometry_ resol_;
  const Variables ctlvars_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostFct3DVar<MODEL, OBS>::CostFct3DVar(const Parameters_ & params,
                                       const eckit::mpi::Comm & comm)
  : CostFunction<MODEL, OBS>::CostFunction(),
    windowLength_(), windowHalf_(), comm_(comm),
    resol_(params.geometry, comm),
    ctlvars_(params.analysisVariables)
{
  Log::trace() << "CostFct3DVar::CostFct3DVar start" << std::endl;
  windowLength_ = params.windowLength;
  windowBegin_ = params.windowBegin;
  windowEnd_ = windowBegin_ + windowLength_;
  windowHalf_ = windowBegin_ + windowLength_/2;

  this->setupTerms(params.toConfiguration());  // Background is read here

  Log::info() << "3DVar window: begin = " << windowBegin_ << ", end = " << windowEnd_ << std::endl;
  Log::trace() << "CostFct3DVar::CostFct3DVar done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJb3D<MODEL> * CostFct3DVar<MODEL, OBS>::newJb(const eckit::Configuration & jbConf,
                                                  const Geometry_ & resol,
                                                  const CtrlVar_ & xb) const {
  Log::trace() << "CostFct3DVar::newJb" << std::endl;
  return new CostJb3D<MODEL>(jbConf, resol, ctlvars_, util::Duration(0), xb.state());
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJo<MODEL, OBS> * CostFct3DVar<MODEL, OBS>::newJo(
    const ObserversParameters<MODEL, OBS> & joParams) const {
  Log::trace() << "CostFct3DVar::newJo" << std::endl;
  return new CostJo<MODEL, OBS>(joParams, comm_, windowBegin_, windowEnd_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostTermBase<MODEL, OBS> * CostFct3DVar<MODEL, OBS>::newJc(const eckit::Configuration & jcConf,
                                                           const Geometry_ &) const {
  Log::trace() << "CostFct3DVar::newJc" << std::endl;
// For now there is no Jc that can work with 3D-Var
  CostTermBase<MODEL, OBS> * pjc = 0;
  return pjc;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct3DVar<MODEL, OBS>::runNL(CtrlVar_ & xx, PostProcessor<State_> & post) const {
  Log::trace() << "CostFct3DVar::runNL start" << std::endl;
  ASSERT(xx.state().validTime() == windowHalf_);

  post.initialize(xx.state(), windowHalf_, windowLength_);
  post.process(xx.state());
  post.finalize(xx.state());

  ASSERT(xx.state().validTime() == windowHalf_);
  Log::trace() << "CostFct3DVar::runNL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFct3DVar<MODEL, OBS>::doLinearize(const Geometry_ & res, const eckit::Configuration & conf,
                                           const CtrlVar_ &, const CtrlVar_ &,
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
  ASSERT(dx.state().validTime() == windowHalf_);

  cost.initializeTL(dx.state(), windowHalf_, windowLength_);
  post.initialize(dx.state(), windowHalf_, windowLength_);

  cost.processTL(dx.state());
  post.process(dx.state());

  cost.finalizeTL(dx.state());
  post.finalize(dx.state());

  ASSERT(dx.state().validTime() == windowHalf_);
  Log::trace() << "CostFct3DVar::runTLM done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct3DVar<MODEL, OBS>::zeroAD(CtrlInc_ & dx) const {
  Log::trace() << "CostFct3DVar::zeroAD start" << std::endl;
  dx.state().zero(windowHalf_);
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
  ASSERT(dx.state().validTime() == windowHalf_);

  post.initialize(dx.state(), windowHalf_, windowLength_);
  cost.initializeAD(dx.state(), windowHalf_, windowLength_);

  cost.processAD(dx.state());
  post.process(dx.state());

  cost.finalizeAD(dx.state());
  post.finalize(dx.state());

  ASSERT(dx.state().validTime() == windowHalf_);

  Log::trace() << "CostFct3DVar::runADJ done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFct3DVar<MODEL, OBS>::addIncr(CtrlVar_ & xx, const CtrlInc_ & dx,
                                       PostProcessor<Increment_> &) const {
  Log::trace() << "CostFct3DVar::addIncr start" << std::endl;
  ASSERT(xx.state().validTime() == windowHalf_);
  ASSERT(dx.state().validTime() == windowHalf_);
  xx.state() += dx.state();
  Log::trace() << "CostFct3DVar::addIncr done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFCT3DVAR_H_
