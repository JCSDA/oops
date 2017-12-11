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

#include "eckit/config/Configuration.h"
#include "util/Logger.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJb3D.h"
#include "oops/assimilation/CostJo.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTL.h"
#include "oops/base/PostProcessorAD.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Model.h"
#include "oops/interface/State.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace oops {

/// 3D-Var Cost Function
/*!
 * This class is not really necessary since it is only a special
 * case of the more general 4D-Var cost function. It is provided
 * for readability.
 */

// -----------------------------------------------------------------------------

template<typename MODEL> class CostFct3DVar : public CostFunction<MODEL> {
  typedef Increment<MODEL>           Increment_;
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef ControlVariable<MODEL>     CtrlVar_;
  typedef CostFunction<MODEL>        CostFct_;
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;
  typedef Model<MODEL>               Model_;

 public:
  CostFct3DVar(const eckit::Configuration &, const Geometry_ &, const Model_ &);
  ~CostFct3DVar() {}

  void runTLM(CtrlInc_ &, PostProcessorTL<Increment_> &,
              PostProcessor<Increment_>, 
              const bool idModel = false) const override;
  void runADJ(CtrlInc_ &, PostProcessorAD<Increment_> &,
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

  util::Duration windowLength_;
  util::DateTime windowBegin_;
  util::DateTime windowEnd_;
  util::DateTime windowHalf_;
  util::Duration zero_;
  const Variables ctlvars_;
};

// =============================================================================

template<typename MODEL>
CostFct3DVar<MODEL>::CostFct3DVar(const eckit::Configuration & config,
                                  const Geometry_ & resol, const Model_ & model)
  : CostFunction<MODEL>::CostFunction(resol, model),
    windowLength_(), windowHalf_(), zero_(0), ctlvars_(config)
{
  windowLength_ = util::Duration(config.getString("window_length"));
  windowBegin_ = util::DateTime(config.getString("window_begin"));
  windowEnd_ = windowBegin_ + windowLength_;
  windowHalf_ = windowBegin_ + windowLength_/2;
  this->setupTerms(config);
  Log::trace() << "CostFct3DVar constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
CostJb3D<MODEL> * CostFct3DVar<MODEL>::newJb(const eckit::Configuration & jbConf,
                                             const Geometry_ & resol,
                                             const CtrlVar_ & xb) const {
  ASSERT(xb.state().checkStatesNumber(1));
  return new CostJb3D<MODEL>(jbConf, resol, ctlvars_, zero_, xb.state()[0]);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
CostJo<MODEL> * CostFct3DVar<MODEL>::newJo(const eckit::Configuration & joConf) const {
  return new CostJo<MODEL>(joConf, windowBegin_, windowEnd_, windowLength_);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
CostTermBase<MODEL> * CostFct3DVar<MODEL>::newJc(const eckit::Configuration & jcConf,
                                                 const Geometry_ &) const {
// For now there is no Jc that can work with 3D-Var
  CostTermBase<MODEL> * pjc = 0;
  return pjc;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct3DVar<MODEL>::runNL(CtrlVar_ & xx,
                                PostProcessor<State_> & post) const {
  ASSERT(xx.state().checkStatesNumber(1));
  ASSERT(xx.state()[0].validTime() == windowHalf_);
  CostFct_::getModel().forecast(xx.state()[0], xx.modVar(), util::Duration(0), post);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct3DVar<MODEL>::runTLM(CtrlInc_ & dx,
                                 PostProcessorTL<Increment_> & cost,
                                 PostProcessor<Increment_> post,
                                 const bool idModel) const {
  ASSERT(dx.state()[0].validTime() == windowHalf_);
  CostFct_::getTLM().forecastTL(dx.state()[0], dx.modVar(), util::Duration(0), post, cost);
  ASSERT(dx.state()[0].validTime() == windowHalf_);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct3DVar<MODEL>::zeroAD(CtrlInc_ & dx) const {
  dx.state()[0].zero(windowHalf_);
  dx.modVar().zero();
  dx.obsVar().zero();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct3DVar<MODEL>::runADJ(CtrlInc_ & dx,
                                 PostProcessorAD<Increment_> & cost,
                                 PostProcessor<Increment_> post,
                                 const bool idModel) const {
  ASSERT(dx.state()[0].validTime() == windowHalf_);
  CostFct_::getTLM().forecastAD(dx.state()[0], dx.modVar(), util::Duration(0), post, cost);
  ASSERT(dx.state()[0].validTime() == windowHalf_);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostFct3DVar<MODEL>::addIncr(CtrlVar_ & xx, const CtrlInc_ & dx,
                                  PostProcessor<Increment_> &) const {
  ASSERT(xx.state().checkStatesNumber(1));
  ASSERT(xx.state()[0].validTime() == windowHalf_);
  ASSERT(dx.state()[0].validTime() == windowHalf_);
  xx.state()[0] += dx.state()[0];
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFCT3DVAR_H_
