/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_COSTFCT4DENSVAR_H_
#define OOPS_ASSIMILATION_COSTFCT4DENSVAR_H_

#include <map>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJb4D.h"
#include "oops/assimilation/CostJcDFI.h"
#include "oops/assimilation/CostJo.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/StateInfo.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Model.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

/// 4D-Ens-Var Cost Function
/*!
 *  Although so far only used for 4D-Ens-Var this cost function can
 *  be interpreted more generally as a four dimensional 3D-Var in the
 *  sense that the control variable is 4D (like weak-constraint 4D-Var)
 *  but the observation operator is 3D (does not involve the forecast
 *  model).
 */

// -----------------------------------------------------------------------------

template<typename MODEL> class CostFct4DEnsVar : public CostFunction<MODEL> {
  typedef Increment<MODEL>           Increment_;
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef ControlVariable<MODEL>     CtrlVar_;
  typedef CostFunction<MODEL>        CostFct_;
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;
  typedef Model<MODEL>               Model_;

 public:
  CostFct4DEnsVar(const eckit::Configuration &, const Geometry_ &, const Model_ &);
  ~CostFct4DEnsVar() {}

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

  CostJb4D<MODEL>     * newJb(const eckit::Configuration &, const Geometry_ &,
                              const CtrlVar_ &) const override;
  CostJo<MODEL>       * newJo(const eckit::Configuration &) const override;
  CostTermBase<MODEL> * newJc(const eckit::Configuration &, const Geometry_ &) const override;

  util::Duration windowLength_;
  util::DateTime windowBegin_;
  util::DateTime windowEnd_;
  util::Duration windowSub_;
  util::Duration zero_;
  unsigned int ncontrol_;
  const Variables ctlvars_;
};

// =============================================================================

template<typename MODEL>
CostFct4DEnsVar<MODEL>::CostFct4DEnsVar(const eckit::Configuration & config,
                                        const Geometry_ & resol, const Model_ & model)
  : CostFunction<MODEL>::CostFunction(config, resol, model), zero_(0), ctlvars_(config)
{
  windowLength_ = util::Duration(config.getString("window_length"));
  windowBegin_ = util::DateTime(config.getString("window_begin"));
  windowEnd_ = windowBegin_ + windowLength_;
  windowSub_ = util::Duration(config.getString("window_sub"));

  ncontrol_ = windowLength_.toSeconds() / windowSub_.toSeconds();
  ASSERT(windowLength_.toSeconds() == windowSub_.toSeconds()*ncontrol_);

  this->setupTerms(config);
  Log::trace() << "CostFct4DEnsVar constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
CostJb4D<MODEL> * CostFct4DEnsVar<MODEL>::newJb(const eckit::Configuration & jbConf,
                                                const Geometry_ & resol,
                                                const CtrlVar_ & xb) const {
  return new CostJb4D<MODEL>(jbConf, resol, ctlvars_, zero_, xb.state());
}

// -----------------------------------------------------------------------------

template <typename MODEL>
CostJo<MODEL> * CostFct4DEnsVar<MODEL>::newJo(const eckit::Configuration & joConf) const {
  return new CostJo<MODEL>(joConf, windowBegin_, windowEnd_, windowSub_, true);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
CostTermBase<MODEL> * CostFct4DEnsVar<MODEL>::newJc(const eckit::Configuration & jcConf,
                                                    const Geometry_ & resol) const {
  const eckit::LocalConfiguration jcdfi(jcConf, "jcdfi");
  const util::DateTime vt(windowBegin_ + windowLength_/2);
  return new CostJcDFI<MODEL>(jcdfi, resol, vt, windowLength_, windowSub_);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct4DEnsVar<MODEL>::runNL(CtrlVar_ & xx,
                                   PostProcessor<State_> & post) const {
  State_ xm(xx.state()[0].geometry(), CostFct_::getModel().variables(), windowBegin_);
  for (unsigned int jsub = 0; jsub <= ncontrol_; ++jsub) {
    util::DateTime now(windowBegin_ + jsub*windowSub_);

    ASSERT(xx.state()[jsub].validTime() == now);
    xm = xx.state()[jsub];
    CostFct_::getModel().forecast(xm, xx.modVar(), util::Duration(0), post);
    xx.state()[jsub] = xm;
    ASSERT(xx.state()[jsub].validTime() == now);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct4DEnsVar<MODEL>::runTLM(CtrlInc_ & dx,
                                    PostProcessorTLAD<MODEL> & cost,
                                    PostProcessor<Increment_> post,
                                    const bool idModel) const {
  for (unsigned int jsub = 0; jsub <= ncontrol_; ++jsub) {
    util::DateTime now(windowBegin_ + jsub*windowSub_);

    ASSERT(dx.state()[jsub].validTime() == now);
    CostFct_::getTLM(jsub).forecastTL(dx.state()[jsub], dx.modVar(), util::Duration(0), post, cost);
    ASSERT(dx.state()[jsub].validTime() == now);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct4DEnsVar<MODEL>::zeroAD(CtrlInc_ & dx) const {
  util::DateTime now(windowBegin_);
  for (unsigned int jsub = 0; jsub <= ncontrol_; ++jsub) {
    dx.state()[jsub].zero(now);
    now += windowSub_;
  }
  dx.modVar().zero();
  dx.obsVar().zero();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFct4DEnsVar<MODEL>::runADJ(CtrlInc_ & dx,
                                    PostProcessorTLAD<MODEL> & cost,
                                    PostProcessor<Increment_> post,
                                    const bool idModel) const {
  for (int jsub = ncontrol_; jsub >= 0; --jsub) {
    util::DateTime now(windowBegin_ + jsub*windowSub_);

    ASSERT(dx.state()[jsub].validTime() == now);
    CostFct_::getTLM(jsub).forecastAD(dx.state()[jsub], dx.modVar(), util::Duration(0), post, cost);
    ASSERT(dx.state()[jsub].validTime() == now);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostFct4DEnsVar<MODEL>::addIncr(CtrlVar_ & xx, const CtrlInc_ & dx,
                                     PostProcessor<Increment_> &) const {
  for (unsigned jsub = 0; jsub <= ncontrol_; ++jsub) {
    xx.state()[jsub] += dx.state()[jsub];
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFCT4DENSVAR_H_
