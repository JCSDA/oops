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
#include "eckit/mpi/Comm.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJbJq.h"
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

/// Weak Constraint 4D-Var Cost Function
/*!
 * General weak constraint constraint 4D-Var cost function.
 */

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS> class CostFctWeak : public CostFunction<MODEL, OBS> {
  typedef Increment<MODEL>                Increment_;
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef ControlVariable<MODEL, OBS>     CtrlVar_;
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef Geometry<MODEL>                 Geometry_;
  typedef State<MODEL>                    State_;
  typedef Model<MODEL>                    Model_;
  typedef VariableChangeBase<MODEL>       VarCha_;
  typedef LinearVariableChangeBase<MODEL> LinVarCha_;

 public:
  CostFctWeak(const eckit::Configuration &, const eckit::mpi::Comm &);
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
  CostJo<MODEL, OBS>       * newJo(const eckit::Configuration &) const override;
  CostTermBase<MODEL, OBS> * newJc(const eckit::Configuration &, const Geometry_ &) const override;
  void doLinearize(const Geometry_ &, const eckit::Configuration &,
                   const CtrlVar_ &, const CtrlVar_ &) override;
  const Geometry_ & geometry() const override {return resol_;}

  util::Duration windowLength_;
  util::DateTime windowBegin_;
  util::DateTime windowEnd_;
  util::Duration windowSub_;
  unsigned int nsubwin_;
  const eckit::mpi::Comm & comm_;
  Geometry_ resol_;
  Model_ model_;
  bool tlforcing_;
  const Variables ctlvars_;
  std::unique_ptr<VarCha_> an2model_;
  std::unique_ptr<LinVarCha_> inc2model_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostFctWeak<MODEL, OBS>::CostFctWeak(const eckit::Configuration & config,
                                const eckit::mpi::Comm & comm)
  : CostFunction<MODEL, OBS>::CostFunction(config), comm_(comm),
    resol_(eckit::LocalConfiguration(config, "resolution"), comm),
    model_(resol_, eckit::LocalConfiguration(config, "model")),
    tlforcing_(false), ctlvars_(config), an2model_(), inc2model_()
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

  an2model_.reset(VariableChangeFactory<MODEL>::create(config, resol_));

  Log::trace() << "CostFctWeak constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJbJq<MODEL> * CostFctWeak<MODEL, OBS>::newJb(const eckit::Configuration & jbConf,
                                            const Geometry_ & resol,
                                            const CtrlVar_ & xb) const {
  return new CostJbJq<MODEL>(jbConf, resol, ctlvars_, windowSub_, xb.state(), tlforcing_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJo<MODEL, OBS> * CostFctWeak<MODEL, OBS>::newJo(const eckit::Configuration & joConf) const {
  return new CostJo<MODEL, OBS>(joConf, comm_, windowBegin_, windowEnd_, util::Duration(0), true);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostTermBase<MODEL, OBS> * CostFctWeak<MODEL, OBS>::newJc(const eckit::Configuration & jcConf,
                                                const Geometry_ & resol) const {
  const eckit::LocalConfiguration jcdfi(jcConf, "jcdfi");
  const util::DateTime vt(windowBegin_ + windowLength_/2);
  return new CostJcDFI<MODEL, OBS>(jcdfi, resol, vt, windowLength_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::runNL(CtrlVar_ & xx,
                               PostProcessor<State_> & post) const {
  State_ xm(xx.state()[0].geometry(), model_.variables(), windowBegin_);
  for (unsigned int jsub = 0; jsub < nsubwin_; ++jsub) {
    util::DateTime bgn(windowBegin_ + jsub*windowSub_);
    util::DateTime end(bgn + windowSub_);

    ASSERT(xx.state()[jsub].validTime() == bgn);
    an2model_->changeVar(xx.state()[jsub], xm);
    model_.forecast(xm, xx.modVar(), windowSub_, post);
    an2model_->changeVarInverse(xm, xx.state()[jsub]);
    ASSERT(xx.state()[jsub].validTime() == end);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::doLinearize(const Geometry_ & resol,
                                     const eckit::Configuration & innerConf,
                                     const CtrlVar_ & bg, const CtrlVar_ & fg) {
  Log::trace() << "CostFctWeak::doLinearize start" << std::endl;
  eckit::LocalConfiguration conf(innerConf, "linearmodel");
  inc2model_.reset(LinearVariableChangeFactory<MODEL>::create(bg.state()[0], fg.state()[0],
                                                             resol, conf));
  inc2model_->setInputVariables(ctlvars_);
  inc2model_->setOutputVariables(CostFct_::getTLM().variables());
  Log::trace() << "CostFctWeak::doLinearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::runTLM(CtrlInc_ & dx,
                                PostProcessorTLAD<MODEL> & cost,
                                PostProcessor<Increment_> post,
                                const bool idModel) const {
  Increment_ dxmodel(dx.state()[0].geometry(), CostFct_::getTLM().variables(), windowBegin_);

  for (int jsub = dx.state().first(); jsub <= dx.state().last(); ++jsub) {
    util::DateTime bgn(windowBegin_ + jsub*windowSub_);
    util::DateTime end(bgn + windowSub_);

    ASSERT(dx.state()[jsub].validTime() == bgn);
    if (tlforcing_ && jsub > 0) dx.state()[jsub] += dx.state()[jsub-1];
    inc2model_->multiply(dx.state()[jsub], dxmodel);
    CostFct_::getTLM(jsub).forecastTL(dxmodel, dx.modVar(), windowSub_, post, cost, idModel);
    inc2model_->multiplyInverse(dxmodel, dx.state()[jsub]);
    ASSERT(dx.state()[jsub].validTime() == end);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::runTLM(CtrlInc_ & dx, const bool idModel) const {
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
      inc2model_->multiply(dx.state()[jsub], dxmodel);
      CostFct_::getTLM(jsub).forecastTL(dxmodel, dx.modVar(), windowSub_, post, cost);
      inc2model_->multiplyInverse(dxmodel, dx.state()[jsub]);
    }

    ASSERT(dx.state()[jsub].validTime() == end);
  }

  dx.state().shift_forward();
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::zeroAD(CtrlInc_ & dx) const {
  for (int jsub = dx.state().first(); jsub <= dx.state().last(); ++jsub) {
    util::DateTime end(windowBegin_ + (jsub+1)*windowSub_);
    dx.state()[jsub].zero(end);
  }
  dx.modVar().zero();
  dx.obsVar().zero();
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::runADJ(CtrlInc_ & dx,
                                PostProcessorTLAD<MODEL> & cost,
                                PostProcessor<Increment_> post,
                                const bool idModel) const {
  Increment_ dxmodel(dx.state()[0].geometry(), CostFct_::getTLM().variables(), windowEnd_);
  for (int jsub = dx.state().last(); jsub >= dx.state().first(); --jsub) {
    util::DateTime bgn(windowBegin_ + jsub*windowSub_);
    util::DateTime end(bgn + windowSub_);

    ASSERT(dx.state()[jsub].validTime() == end);
    inc2model_->multiplyInverseAD(dx.state()[jsub], dxmodel);
    CostFct_::getTLM(jsub).forecastAD(dxmodel, dx.modVar(), windowSub_, post, cost, idModel);
    inc2model_->multiplyAD(dxmodel, dx.state()[jsub]);
    if (tlforcing_ && jsub > 0) dx.state()[jsub-1] += dx.state()[jsub];
    ASSERT(dx.state()[jsub].validTime() == bgn);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::runADJ(CtrlInc_ & dx, const bool idModel) const {
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
      inc2model_->multiplyInverseAD(dx.state()[jsub], dxmodel);
      CostFct_::getTLM(jsub).forecastAD(dxmodel, dx.modVar(), windowSub_, post, cost);
      inc2model_->multiplyAD(dxmodel, dx.state()[jsub]);
    }

    ASSERT(dx.state()[jsub].validTime() == bgn);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::addIncr(CtrlVar_ & xx, const CtrlInc_ & dx,
                                 PostProcessor<Increment_> & post) const {
  if (tlforcing_) {
    Increment_ xi(dx.state()[0]);
    Increment_ dxmodel(dx.state()[0].geometry(), CostFct_::getTLM().variables(), windowBegin_);
    for (unsigned int jsub = 0; jsub < nsubwin_; ++jsub) {
      if (jsub > 0) xi += dx.state()[jsub];
      xx.state()[jsub] += xi;
      if (jsub < nsubwin_-1) {
        inc2model_->multiply(xi, dxmodel);
        CostFct_::getTLM(jsub).forecastTL(dxmodel, dx.modVar(), windowSub_, post);
        inc2model_->multiplyInverse(dxmodel, xi);
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
