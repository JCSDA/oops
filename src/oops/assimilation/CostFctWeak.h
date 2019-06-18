/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_COSTFCTWEAK_H_
#define OOPS_ASSIMILATION_COSTFCTWEAK_H_

#include <map>
#include <memory>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJbJq.h"
#include "oops/assimilation/CostJcDFI.h"
#include "oops/assimilation/CostJo.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/LinearVariableChangeBase.h"
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

/// Weak Constraint 4D-Var Cost Function
/*!
 * General weak constraint constraint 4D-Var cost function.
 */

// -----------------------------------------------------------------------------

template<typename MODEL> class CostFctWeak : public CostFunction<MODEL> {
  typedef Increment<MODEL>           Increment_;
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef ControlVariable<MODEL>     CtrlVar_;
  typedef CostFunction<MODEL>        CostFct_;
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;
  typedef Model<MODEL>               Model_;
  typedef LinearVariableChangeBase<MODEL> ChangeVar_;

 public:
  CostFctWeak(const eckit::Configuration &, const Geometry_ &, const Model_ &);
  ~CostFctWeak() {}

  void runTLM(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
              PostProcessor<Increment_>,
              const bool idModel = false) const override;
  void runADJ(CtrlInc_ &, PostProcessorTLAD<MODEL> &,
              PostProcessor<Increment_>,
              const bool idModel = false) const override;
  void zeroAD(CtrlInc_ &) const override;

  void runTLM(CtrlInc_ &, const bool idModel = false) const;
  void runADJ(CtrlInc_ &, const bool idModel = false) const;
  void runNL(CtrlVar_ &, PostProcessor<State_> &) const override;

 private:
  void addIncr(CtrlVar_ &, const CtrlInc_ &, PostProcessor<Increment_> &) const override;

  CostJbJq<MODEL>     * newJb(const eckit::Configuration &, const Geometry_ &,
                              const CtrlVar_ &) const override;
  CostJo<MODEL>       * newJo(const eckit::Configuration &) const override;
  CostTermBase<MODEL> * newJc(const eckit::Configuration &, const Geometry_ &) const override;
  void doLinearize(const Geometry_ &, const eckit::Configuration &,
                   const CtrlVar_ &, const CtrlVar_ &) override;

  util::Duration windowLength_;
  util::DateTime windowBegin_;
  util::DateTime windowEnd_;
  util::Duration windowSub_;
  unsigned int nsubwin_;
  bool tlforcing_;
  const Variables ctlvars_;
  std::unique_ptr<ChangeVar_> an2model_;
};

// =============================================================================

template<typename MODEL>
CostFctWeak<MODEL>::CostFctWeak(const eckit::Configuration & config,
                                const Geometry_ & resol, const Model_ & model)
  : CostFunction<MODEL>::CostFunction(config, resol, model),
    tlforcing_(false), ctlvars_(config), an2model_()
{
  windowLength_ = util::Duration(config.getString("window_length"));
  windowBegin_ = util::DateTime(config.getString("window_begin"));
  windowEnd_ = windowBegin_ + windowLength_;
  windowSub_ = util::Duration(config.getString("window_sub"));

  nsubwin_ = windowLength_.toSeconds() / windowSub_.toSeconds();
  ASSERT(windowLength_.toSeconds() == windowSub_.toSeconds()*nsubwin_);

  if (config.has("tlforcing")) {
    if (config.getString("tlforcing") == "on") tlforcing_ = true;
  }

  this->setupTerms(config);
  Log::trace() << "CostFctWeak constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
CostJbJq<MODEL> * CostFctWeak<MODEL>::newJb(const eckit::Configuration & jbConf,
                                            const Geometry_ & resol,
                                            const CtrlVar_ & xb) const {
  return new CostJbJq<MODEL>(jbConf, resol, ctlvars_, windowSub_, xb.state(), tlforcing_);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
CostJo<MODEL> * CostFctWeak<MODEL>::newJo(const eckit::Configuration & joConf) const {
  return new CostJo<MODEL>(joConf, windowBegin_, windowEnd_, util::Duration(0), true);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
CostTermBase<MODEL> * CostFctWeak<MODEL>::newJc(const eckit::Configuration & jcConf,
                                                const Geometry_ & resol) const {
  const eckit::LocalConfiguration jcdfi(jcConf, "jcdfi");
  const util::DateTime vt(windowBegin_ + windowLength_/2);
  return new CostJcDFI<MODEL>(jcdfi, resol, vt, windowLength_);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFctWeak<MODEL>::runNL(CtrlVar_ & xx,
                               PostProcessor<State_> & post) const {
  State_ xm(xx.state()[0].geometry(), CostFct_::getModel().variables(), windowBegin_);
  for (unsigned int jsub = 0; jsub < nsubwin_; ++jsub) {
    util::DateTime bgn(windowBegin_ + jsub*windowSub_);
    util::DateTime end(bgn + windowSub_);

    ASSERT(xx.state()[jsub].validTime() == bgn);
    xm = xx.state()[jsub];
    CostFct_::getModel().forecast(xm, xx.modVar(), windowSub_, post);
    xx.state()[jsub] = xm;
    ASSERT(xx.state()[jsub].validTime() == end);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostFctWeak<MODEL>::doLinearize(const Geometry_ & resol,
                                     const eckit::Configuration & innerConf,
                                     const CtrlVar_ & bg, const CtrlVar_ & fg) {
  Log::trace() << "CostFctWeak::doLinearize start" << std::endl;
  eckit::LocalConfiguration conf(innerConf, "linearmodel");
  an2model_.reset(LinearVariableChangeFactory<MODEL>::create(bg.state()[0], fg.state()[0],
                                                             resol, conf));
  an2model_->setInputVariables(ctlvars_);
  an2model_->setOutputVariables(CostFct_::getTLM().variables());
  Log::trace() << "CostFctWeak::doLinearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFctWeak<MODEL>::runTLM(CtrlInc_ & dx,
                                PostProcessorTLAD<MODEL> & cost,
                                PostProcessor<Increment_> post,
                                const bool idModel) const {
  Increment_ dxmodel(dx.state()[0].geometry(), CostFct_::getTLM().variables(), windowBegin_);

  for (int jsub = dx.state().first(); jsub <= dx.state().last(); ++jsub) {
    util::DateTime bgn(windowBegin_ + jsub*windowSub_);
    util::DateTime end(bgn + windowSub_);

    ASSERT(dx.state()[jsub].validTime() == bgn);
    if (tlforcing_ && jsub > 0) dx.state()[jsub] += dx.state()[jsub-1];
    an2model_->multiply(dx.state()[jsub], dxmodel);
    CostFct_::getTLM(jsub).forecastTL(dxmodel, dx.modVar(), windowSub_, post, cost, idModel);
    an2model_->multiplyInverse(dxmodel, dx.state()[jsub]);
    ASSERT(dx.state()[jsub].validTime() == end);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFctWeak<MODEL>::runTLM(CtrlInc_ & dx, const bool idModel) const {
  PostProcessor<Increment_> post;
  PostProcessorTLAD<MODEL> cost;
  ASSERT(!tlforcing_);
  Increment_ dxmodel(dx.state()[0].geometry(), CostFct_::getTLM().variables(), windowBegin_);

  for (int jsub = dx.state().first(); jsub <= dx.state().last(); ++jsub) {
    util::DateTime bgn(windowBegin_ + jsub*windowSub_);
    util::DateTime end(bgn + windowSub_);
    ASSERT(dx.state()[jsub].validTime() == bgn);

    if (idModel) {
      dx.state()[jsub].updateTime(windowSub_);
    } else {
      an2model_->multiply(dx.state()[jsub], dxmodel);
      CostFct_::getTLM(jsub).forecastTL(dxmodel, dx.modVar(), windowSub_, post, cost);
      an2model_->multiplyInverse(dxmodel, dx.state()[jsub]);
    }

    ASSERT(dx.state()[jsub].validTime() == end);
  }

  dx.state().shift_forward();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFctWeak<MODEL>::zeroAD(CtrlInc_ & dx) const {
  for (int jsub = dx.state().first(); jsub <= dx.state().last(); ++jsub) {
    util::DateTime end(windowBegin_ + (jsub+1)*windowSub_);
    dx.state()[jsub].zero(end);
  }
  dx.modVar().zero();
  dx.obsVar().zero();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFctWeak<MODEL>::runADJ(CtrlInc_ & dx,
                                PostProcessorTLAD<MODEL> & cost,
                                PostProcessor<Increment_> post,
                                const bool idModel) const {
  Increment_ dxmodel(dx.state()[0].geometry(), CostFct_::getTLM().variables(), windowEnd_);
  for (int jsub = dx.state().last(); jsub >= dx.state().first(); --jsub) {
    util::DateTime bgn(windowBegin_ + jsub*windowSub_);
    util::DateTime end(bgn + windowSub_);

    ASSERT(dx.state()[jsub].validTime() == end);
    an2model_->multiplyInverseAD(dx.state()[jsub], dxmodel);
    CostFct_::getTLM(jsub).forecastAD(dxmodel, dx.modVar(), windowSub_, post, cost, idModel);
    an2model_->multiplyAD(dxmodel, dx.state()[jsub]);
    if (tlforcing_ && jsub > 0) dx.state()[jsub-1] += dx.state()[jsub];
    ASSERT(dx.state()[jsub].validTime() == bgn);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void CostFctWeak<MODEL>::runADJ(CtrlInc_ & dx, const bool idModel) const {
  PostProcessor<Increment_> post;
  PostProcessorTLAD<MODEL> cost;
  ASSERT(!tlforcing_);
  Increment_ dxmodel(dx.state()[0].geometry(), CostFct_::getTLM().variables(), windowEnd_);

  dx.state().shift_backward();

  for (int jsub = dx.state().last(); jsub >= dx.state().first(); --jsub) {
    util::DateTime bgn(windowBegin_ + jsub*windowSub_);
    util::DateTime end(bgn + windowSub_);
    ASSERT(dx.state()[jsub].validTime() == end);

    if (idModel) {
      dx.state()[jsub].updateTime(-windowSub_);
    } else {
      an2model_->multiplyInverseAD(dx.state()[jsub], dxmodel);
      CostFct_::getTLM(jsub).forecastAD(dxmodel, dx.modVar(), windowSub_, post, cost);
      an2model_->multiplyAD(dxmodel, dx.state()[jsub]);
    }

    ASSERT(dx.state()[jsub].validTime() == bgn);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostFctWeak<MODEL>::addIncr(CtrlVar_ & xx, const CtrlInc_ & dx,
                                 PostProcessor<Increment_> & post) const {
  if (tlforcing_) {
    Increment_ xi(dx.state()[0]);
    Increment_ dxmodel(dx.state()[0].geometry(), CostFct_::getTLM().variables(), windowBegin_);
    for (unsigned int jsub = 0; jsub < nsubwin_; ++jsub) {
      if (jsub > 0) xi += dx.state()[jsub];
      xx.state()[jsub] += xi;
      if (jsub < nsubwin_-1) {
        an2model_->multiply(xi, dxmodel);
        CostFct_::getTLM(jsub).forecastTL(dxmodel, dx.modVar(), windowSub_, post);
        an2model_->multiplyInverse(dxmodel, xi);
      }
    }
  } else {
    for (unsigned int jsub = 0; jsub < nsubwin_; ++jsub) {
      xx.state()[jsub] += dx.state()[jsub];
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFCTWEAK_H_
