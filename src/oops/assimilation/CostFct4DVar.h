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

#include "eckit/config/LocalConfiguration.h"
#include "util/Logger.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJb3D.h"
#include "oops/assimilation/CostJcDFI.h"
#include "oops/assimilation/CostJo.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTL.h"
#include "oops/base/PostProcessorAD.h"
#include "oops/base/StateInfo.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Model.h"
#include "oops/interface/State.h"
#include "oops/interface/Variables.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace oops {

/// Strong Constraint 4D-Var Cost Function
/*!
 * This class is not really necessary since it is only a special
 * case of the more general weak constraint 4D-Var cost function
 * with one sub-window. It is provided for readability.
 */

// -----------------------------------------------------------------------------

template<typename MODEL> class CostFct4DVar : public CostFunction<MODEL> {
  typedef Increment<MODEL>           Increment_;
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef ControlVariable<MODEL>     CtrlVar_;
  typedef CostFunction<MODEL>        CostFct_;
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;
  typedef Model<MODEL>               Model_;
  typedef Variables<MODEL>           Variables_;

 public:
  CostFct4DVar(const eckit::Configuration &, const Geometry_ &, const Model_ &);
  ~CostFct4DVar() {}

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
  const Variables_ ctlvars_;
};

// =============================================================================

template<typename MODEL>
CostFct4DVar<MODEL>::CostFct4DVar(const eckit::Configuration & config,
                                  const Geometry_ & resol, const Model_ & model)
  : CostFunction<MODEL>::CostFunction(resol, model), ctlvars_(config)
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
  CostFct_::getModel().forecast(xx.state()[0], xx.modVar(), windowLength_, post);
  ASSERT(xx.state()[0].validTime() == windowEnd_);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct4DVar<MODEL>::runTLM(CtrlInc_ & dx,
                                 PostProcessorTL<Increment_> & cost,
                                 PostProcessor<Increment_> post,
                                 const bool idModel) const {
  ASSERT(dx.state()[0].validTime() == windowBegin_);
  CostFct_::getTLM().forecastTL(dx.state()[0], dx.modVar(), windowLength_, 
                                post, cost, idModel);
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
                                 PostProcessorAD<Increment_> & cost,
                                 PostProcessor<Increment_> post,
                                 const bool idModel) const {
  ASSERT(dx.state()[0].validTime() == windowEnd_);
  CostFct_::getTLM().forecastAD(dx.state()[0], dx.modVar(), windowLength_, 
                                post, cost, idModel);
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
