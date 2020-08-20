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
#include <memory>


#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJb4D.h"
#include "oops/assimilation/CostJcDFI.h"
#include "oops/assimilation/CostJo.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/StateInfo.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
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

template<typename MODEL, typename OBS> class CostFct4DEnsVar : public CostFunction<MODEL, OBS> {
  typedef Increment<MODEL>                Increment_;
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef ControlVariable<MODEL, OBS>     CtrlVar_;
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef Geometry<MODEL>                 Geometry_;
  typedef State<MODEL>                    State_;

 public:
  CostFct4DEnsVar(const eckit::Configuration &, const eckit::mpi::Comm &);
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
  CostJo<MODEL, OBS>       * newJo(const eckit::Configuration &) const override;
  CostTermBase<MODEL, OBS> * newJc(const eckit::Configuration &, const Geometry_ &) const override;
  void doLinearize(const Geometry_ &, const eckit::Configuration &,
                   const CtrlVar_ &, const CtrlVar_ &,
                   PostProcessor<State_> &, PostProcessorTLAD<MODEL> &) override;
  const Geometry_ & geometry() const override {return resol_;}

  util::Duration windowLength_;
  util::DateTime windowBegin_;
  util::DateTime windowEnd_;
  util::Duration windowSub_;
  const eckit::mpi::Comm & comm_;
  Geometry_ resol_;
  unsigned int ncontrol_;
  const Variables ctlvars_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostFct4DEnsVar<MODEL, OBS>::CostFct4DEnsVar(const eckit::Configuration & config,
                                             const eckit::mpi::Comm & comm)
  : CostFunction<MODEL, OBS>::CostFunction(config), comm_(comm),
    resol_(eckit::LocalConfiguration(config, "geometry"), comm),
    ctlvars_(config, "analysis variables")
{
  windowLength_ = util::Duration(config.getString("window length"));
  windowBegin_ = util::DateTime(config.getString("window begin"));
  windowEnd_ = windowBegin_ + windowLength_;
  windowSub_ = util::Duration(config.getString("subwindow"));

  ncontrol_ = windowLength_.toSeconds() / windowSub_.toSeconds();
  ASSERT(windowLength_.toSeconds() == windowSub_.toSeconds()*ncontrol_);

  this->setupTerms(config);

  Log::trace() << "CostFct4DEnsVar constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJb4D<MODEL> * CostFct4DEnsVar<MODEL, OBS>::newJb(const eckit::Configuration & jbConf,
                                                     const Geometry_ & resol,
                                                     const CtrlVar_ & xb) const {
  return new CostJb4D<MODEL>(jbConf, resol, ctlvars_, xb.state());
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJo<MODEL, OBS> * CostFct4DEnsVar<MODEL, OBS>::newJo(const eckit::Configuration & joConf) const {
  return new CostJo<MODEL, OBS>(joConf, comm_, windowBegin_, windowEnd_, windowSub_, true);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostTermBase<MODEL, OBS> * CostFct4DEnsVar<MODEL, OBS>::newJc(const eckit::Configuration & jcConf,
                                                              const Geometry_ & resol) const {
  const eckit::LocalConfiguration jcdfi(jcConf, "jcdfi");
  const util::DateTime vt(windowBegin_ + windowLength_/2);
  return new CostJcDFI<MODEL, OBS>(jcdfi, resol, vt, windowLength_, windowSub_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct4DEnsVar<MODEL, OBS>::runNL(CtrlVar_ & xx, PostProcessor<State_> & post) const {
  post.initialize(xx.state()[0], windowEnd_, windowSub_);
  for (unsigned int jsub = 0; jsub <= ncontrol_; ++jsub) {
    util::DateTime now(windowBegin_ + jsub*windowSub_);
    ASSERT(xx.state()[jsub].validTime() == now);
    post.process(xx.state()[jsub]);
  }
  post.finalize(xx.state()[ncontrol_]);
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFct4DEnsVar<MODEL, OBS>::doLinearize(const Geometry_ & resol,
                                              const eckit::Configuration & conf,
                                              const CtrlVar_ &, const CtrlVar_ &,
                                              PostProcessor<State_> & pp,
                                              PostProcessorTLAD<MODEL> & pptraj) {
  Log::trace() << "CostFct4DEnsVar::doLinearize start" << std::endl;
  pp.enrollProcessor(new TrajectorySaver<MODEL>(conf, resol, pptraj));
  Log::trace() << "CostFct4DEnsVar::doLinearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct4DEnsVar<MODEL, OBS>::runTLM(CtrlInc_ & dx,
                                         PostProcessorTLAD<MODEL> & cost,
                                         PostProcessor<Increment_> post,
                                         const bool) const {
  cost.initializeTL(dx.state()[0], windowEnd_, windowSub_);
  post.initialize(dx.state()[0], windowEnd_, windowSub_);

  for (unsigned int jsub = 0; jsub <= ncontrol_; ++jsub) {
    util::DateTime now(windowBegin_ + jsub*windowSub_);
    ASSERT(dx.state()[jsub].validTime() == now);

    cost.processTL(dx.state()[jsub]);
    post.process(dx.state()[jsub]);
  }
  cost.finalizeTL(dx.state()[ncontrol_]);
  post.finalize(dx.state()[ncontrol_]);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct4DEnsVar<MODEL, OBS>::zeroAD(CtrlInc_ & dx) const {
  util::DateTime now(windowBegin_);
  for (unsigned int jsub = 0; jsub <= ncontrol_; ++jsub) {
    dx.state()[jsub].zero(now);
    now += windowSub_;
  }
  dx.modVar().zero();
  dx.obsVar().zero();
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct4DEnsVar<MODEL, OBS>::runADJ(CtrlInc_ & dx,
                                         PostProcessorTLAD<MODEL> & cost,
                                         PostProcessor<Increment_> post,
                                         const bool) const {
  post.initialize(dx.state()[ncontrol_], windowBegin_, windowSub_);
  cost.initializeAD(dx.state()[ncontrol_], windowBegin_, windowSub_);

  for (int jsub = ncontrol_; jsub >= 0; --jsub) {
    cost.processAD(dx.state()[jsub]);
    post.process(dx.state()[jsub]);

    util::DateTime now(windowBegin_ + jsub*windowSub_);
    ASSERT(dx.state()[jsub].validTime() == now);
  }

  cost.finalizeAD(dx.state()[0]);
  post.finalize(dx.state()[0]);
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFct4DEnsVar<MODEL, OBS>::addIncr(CtrlVar_ & xx, const CtrlInc_ & dx,
                                          PostProcessor<Increment_> &) const {
  for (unsigned jsub = 0; jsub <= ncontrol_; ++jsub) {
    xx.state()[jsub] += dx.state()[jsub];
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFCT4DENSVAR_H_
