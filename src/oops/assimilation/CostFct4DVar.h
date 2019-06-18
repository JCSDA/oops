/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_COSTFCT4DVAR_H_
#define OOPS_ASSIMILATION_COSTFCT4DVAR_H_

#include <memory>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJb3D.h"
#include "oops/assimilation/CostJcDFI.h"
#include "oops/assimilation/CostJo.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/LinearVariableChangeBase.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/StateInfo.h"
#include "oops/base/VariableChangeBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Model.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

/// Strong Constraint 4D-Var Cost Function
/*!
 * This class is not really necessary since it is only a special
 * case of the more general weak constraint 4D-Var cost function
 * with one sub-window. It is provided for readability.
 */

// -----------------------------------------------------------------------------

template<typename MODEL> class CostFct4DVar : public CostFunction<MODEL> {
  typedef Increment<MODEL>                Increment_;
  typedef ControlIncrement<MODEL>         CtrlInc_;
  typedef ControlVariable<MODEL>          CtrlVar_;
  typedef CostFunction<MODEL>             CostFct_;
  typedef Geometry<MODEL>                 Geometry_;
  typedef State<MODEL>                    State_;
  typedef Model<MODEL>                    Model_;
  typedef LinearVariableChangeBase<MODEL> LinVarCha_;
  typedef VariableChangeBase<MODEL>       VarCha_;

 public:
  CostFct4DVar(const eckit::Configuration &, const Geometry_ &, const Model_ &);
  ~CostFct4DVar() {}

  void runTLM(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
              PostProcessor<Increment_>,
              const bool idModel = false) const override;
  void runADJ(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
              PostProcessor<Increment_>,
              const bool idModel = false) const override;
  void zeroAD(CtrlInc_ &) const override;

  void runNL(CtrlVar_ &, PostProcessor<State_>&) const override;

 private:
  void addIncr(CtrlVar_ &, const CtrlInc_ &, PostProcessor<Increment_>&) const override;

  CostJb3D<MODEL>     * newJb(const eckit::Configuration &, const Geometry_ &,
                              const CtrlVar_ &) const override;
  CostJo<MODEL>       * newJo(const eckit::Configuration &) const override;
  CostTermBase<MODEL> * newJc(const eckit::Configuration &, const Geometry_ &) const override;
  void doLinearize(const Geometry_ &, const eckit::Configuration &,
                   const CtrlVar_ &, const CtrlVar_ &) override;

  util::Duration windowLength_;
  util::DateTime windowBegin_;
  util::DateTime windowEnd_;
  const Variables ctlvars_;
  std::unique_ptr<VarCha_> an2model_;
  std::unique_ptr<LinVarCha_> inc2model_;
};

// =============================================================================

template<typename MODEL>
CostFct4DVar<MODEL>::CostFct4DVar(const eckit::Configuration & config,
                                  const Geometry_ & resol, const Model_ & model)
  : CostFunction<MODEL>::CostFunction(config, resol, model), ctlvars_(config),
    an2model_(VariableChangeFactory<MODEL>::create(config, resol)), inc2model_()
{
  Log::trace() << "CostFct4DVar:CostFct4DVar" << std::endl;
  windowLength_ = util::Duration(config.getString("window_length"));
  windowBegin_ = util::DateTime(config.getString("window_begin"));
  windowEnd_ = windowBegin_ + windowLength_;
  this->setupTerms(config);
  Log::trace() << "CostFct4DVar constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
CostJb3D<MODEL> * CostFct4DVar<MODEL>::newJb(const eckit::Configuration & jbConf,
                                             const Geometry_ & resol,
                                             const CtrlVar_ & xb) const {
  ASSERT(xb.state().checkStatesNumber(1));
  return new CostJb3D<MODEL>(jbConf, resol, ctlvars_, windowLength_, xb.state()[0]);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
CostJo<MODEL> * CostFct4DVar<MODEL>::newJo(const eckit::Configuration & joConf) const {
  return new CostJo<MODEL>(joConf, windowBegin_, windowEnd_, util::Duration(0));
}

// -----------------------------------------------------------------------------

template <typename MODEL>
CostTermBase<MODEL> * CostFct4DVar<MODEL>::newJc(const eckit::Configuration & jcConf,
                                                 const Geometry_ & resol) const {
  const eckit::LocalConfiguration jcdfi(jcConf, "jcdfi");
  const util::DateTime vt(windowBegin_ + windowLength_/2);
  return new CostJcDFI<MODEL>(jcdfi, resol, vt, windowLength_);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct4DVar<MODEL>::runNL(CtrlVar_ & xx,
                                PostProcessor<State_> & post) const {
  ASSERT(xx.state().checkStatesNumber(1));
  ASSERT(xx.state()[0].validTime() == windowBegin_);
  State_ xm(xx.state()[0].geometry(), CostFct_::getModel().variables(), windowBegin_);
  an2model_->changeVar(xx.state()[0], xm);
  CostFct_::getModel().forecast(xm, xx.modVar(), windowLength_, post);
  an2model_->changeVarInverse(xm, xx.state()[0]);
  ASSERT(xx.state()[0].validTime() == windowEnd_);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostFct4DVar<MODEL>::doLinearize(const Geometry_ & resol,
                                      const eckit::Configuration & innerConf,
                                      const CtrlVar_ & bg, const CtrlVar_ & fg) {
  Log::trace() << "CostFct4DVar::doLinearize start" << std::endl;
  eckit::LocalConfiguration conf(innerConf, "linearmodel");
  inc2model_.reset(LinearVariableChangeFactory<MODEL>::create(bg.state()[0], fg.state()[0],
                                                             resol, conf));
  inc2model_->setInputVariables(ctlvars_);
  inc2model_->setOutputVariables(CostFct_::getTLM().variables());
  Log::trace() << "CostFct4DVar::doLinearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct4DVar<MODEL>::runTLM(CtrlInc_ & dx,
                                 PostProcessorTLAD<MODEL> & cost,
                                 PostProcessor<Increment_> post,
                                 const bool idModel) const {
  ASSERT(dx.state()[0].validTime() == windowBegin_);
  Increment_ dxmodel(dx.state()[0].geometry(), CostFct_::getTLM().variables(), windowBegin_);
  inc2model_->multiply(dx.state()[0], dxmodel);
  CostFct_::getTLM().forecastTL(dxmodel, dx.modVar(), windowLength_, post, cost, idModel);
  inc2model_->multiplyInverse(dxmodel, dx.state()[0]);
  ASSERT(dx.state()[0].validTime() == windowEnd_);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct4DVar<MODEL>::zeroAD(CtrlInc_ & dx) const {
  dx.state()[0].zero(windowEnd_);
  dx.modVar().zero();
  dx.obsVar().zero();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct4DVar<MODEL>::runADJ(CtrlInc_ & dx,
                                 PostProcessorTLAD<MODEL> & cost,
                                 PostProcessor<Increment_> post,
                                 const bool idModel) const {
  ASSERT(dx.state()[0].validTime() == windowEnd_);
  Increment_ dxmodel(dx.state()[0].geometry(), CostFct_::getTLM().variables(), windowEnd_);
  inc2model_->multiplyInverseAD(dx.state()[0], dxmodel);
  CostFct_::getTLM().forecastAD(dxmodel, dx.modVar(), windowLength_, post, cost, idModel);
  inc2model_->multiplyAD(dxmodel, dx.state()[0]);
  ASSERT(dx.state()[0].validTime() == windowBegin_);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostFct4DVar<MODEL>::addIncr(CtrlVar_ & xx, const CtrlInc_ & dx,
                                  PostProcessor<Increment_> &) const {
  ASSERT(xx.state().checkStatesNumber(1));
  xx.state()[0] += dx.state()[0];
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFCT4DVAR_H_
