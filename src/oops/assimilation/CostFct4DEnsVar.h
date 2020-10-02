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
#include <string>

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
  const Geometry_ & geometry() const override {return *resol_;}

  util::Duration subWinLength_;
  util::DateTime subWinTime_;
  util::DateTime subWinBegin_;
  util::DateTime subWinEnd_;
  size_t nsubwin_;
  size_t mysubwin_;
  std::unique_ptr<Geometry_> resol_;
  const Variables ctlvars_;
  eckit::mpi::Comm * commSpace_;
  eckit::mpi::Comm * commTime_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostFct4DEnsVar<MODEL, OBS>::CostFct4DEnsVar(const eckit::Configuration & config,
                                             const eckit::mpi::Comm & comm)
  : CostFunction<MODEL, OBS>::CostFunction(config),
    resol_(), ctlvars_(config, "analysis variables")
{
  Log::trace() << "CostFct4DEnsVar::CostFct4DEnsVar start" << std::endl;
  util::Duration windowLength(config.getString("window length"));
  util::DateTime windowBegin(config.getString("window begin"));
  util::DateTime windowEnd = windowBegin + windowLength;
  subWinLength_ = util::Duration(config.getString("subwindow"));

  nsubwin_ = windowLength.toSeconds() / subWinLength_.toSeconds() + 1;  // Not like WC
  ASSERT(windowLength.toSeconds() == subWinLength_.toSeconds() * (nsubwin_ - 1));

  size_t ntasks = comm.size();
  ASSERT(ntasks % nsubwin_ == 0);
  size_t myrank = comm.rank();
  size_t ntaskpslot = ntasks / nsubwin_;
  size_t mysubwin_ = myrank / ntaskpslot;

// Define local sub-window
  subWinTime_  = windowBegin + mysubwin_ * subWinLength_;
  subWinBegin_ = subWinTime_ - subWinLength_/2;
  subWinEnd_   = subWinTime_ + subWinLength_/2;
  if (mysubwin_ == 0) subWinBegin_ = subWinTime_;
  if (mysubwin_ == nsubwin_ - 1) subWinEnd_ = subWinTime_;
  ASSERT(subWinBegin_ >= windowBegin);
  ASSERT(subWinEnd_ <= windowEnd);

// Create a communicator for same sub-window, to be used for communications in space
  std::string sgeom = "comm_geom_" + std::to_string(mysubwin_);
  char const *geomName = sgeom.c_str();
  commSpace_ = &comm.split(mysubwin_, geomName);

// Create a communicator for same local area, to be used for communications in time
  size_t myarea = commSpace_->rank();
  std::string stime = "comm_time_" + std::to_string(myarea);
  char const *timeName = stime.c_str();
  commTime_ = &comm.split(myarea, timeName);
  ASSERT(commTime_->size() == nsubwin_);

// Now can setup the rest
  resol_.reset(new Geometry_(eckit::LocalConfiguration(config, "geometry"),
                             *commSpace_, *commTime_));

  this->setupTerms(config);

  Log::trace() << "CostFct4DEnsVar::CostFct4DEnsVar done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJb4D<MODEL> * CostFct4DEnsVar<MODEL, OBS>::newJb(const eckit::Configuration & jbConf,
                                                     const Geometry_ & resol,
                                                     const CtrlVar_ & xb) const {
  Log::trace() << "CostFct4DEnsVar::newJb" << std::endl;
  return new CostJb4D<MODEL>(jbConf, *commTime_, resol, ctlvars_, xb.state());
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJo<MODEL, OBS> * CostFct4DEnsVar<MODEL, OBS>::newJo(const eckit::Configuration & joConf) const {
  Log::trace() << "CostFct4DEnsVar::newJo" << std::endl;
  return new CostJo<MODEL, OBS>(joConf, *commSpace_,
                                subWinBegin_, subWinEnd_, subWinLength_, *commTime_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostTermBase<MODEL, OBS> * CostFct4DEnsVar<MODEL, OBS>::newJc(const eckit::Configuration & jcConf,
                                                              const Geometry_ & resol) const {
  Log::trace() << "CostFct4DEnsVar::newJc" << std::endl;
//  const eckit::LocalConfiguration jcdfi(jcConf, "jcdfi");
//  const util::DateTime vt(subWinBegin_ + windowLength_/2);
//  return new CostJcDFI<MODEL, OBS>(jcdfi, resol, vt, windowLength_, subWinLength_);
  Log::warning() << "CostFct4DEnsVar::newJc NO Jc" << std::endl;
  CostTermBase<MODEL, OBS> * pjc = 0;
  return pjc;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct4DEnsVar<MODEL, OBS>::runNL(CtrlVar_ & xx, PostProcessor<State_> & post) const {
  Log::trace() << "CostFct4DEnsVar::runNL start" << std::endl;
  ASSERT(xx.state().validTime() == subWinTime_);

  post.initialize(xx.state(), subWinTime_, subWinLength_);
  post.process(xx.state());
  post.finalize(xx.state());

  Log::info() << "CostFct4DEnsVar::runNL: " << xx << std::endl;
  Log::trace() << "CostFct4DEnsVar::runNL done" << std::endl;
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
  Log::trace() << "CostFct4DEnsVar::runTLM start" << std::endl;
  ASSERT(dx.state().validTime() == subWinTime_);

  cost.initializeTL(dx.state(), subWinTime_, subWinLength_);
  post.initialize(dx.state(), subWinTime_, subWinLength_);

  cost.processTL(dx.state());
  post.process(dx.state());

  cost.finalizeTL(dx.state());
  post.finalize(dx.state());

  Log::info() << "CostFct4DEnsVar::runTLM: " << dx << std::endl;
  Log::trace() << "CostFct4DEnsVar::runTLM done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct4DEnsVar<MODEL, OBS>::zeroAD(CtrlInc_ & dx) const {
  Log::trace() << "CostFct4DEnsVar::zeroAD start" << std::endl;
  dx.state().zero(subWinTime_);
  dx.modVar().zero();
  dx.obsVar().zero();
  Log::trace() << "CostFct4DEnsVar::zeroAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct4DEnsVar<MODEL, OBS>::runADJ(CtrlInc_ & dx,
                                         PostProcessorTLAD<MODEL> & cost,
                                         PostProcessor<Increment_> post,
                                         const bool) const {
  Log::trace() << "CostFct4DEnsVar::runADJ start" << std::endl;

  post.initialize(dx.state(), subWinTime_, subWinLength_);
  cost.initializeAD(dx.state(), subWinTime_, subWinLength_);

  cost.processAD(dx.state());
  post.process(dx.state());

  cost.finalizeAD(dx.state());
  post.finalize(dx.state());

  Log::info() << "CostFct4DEnsVar::runADJ: " << dx << std::endl;
  Log::trace() << "CostFct4DEnsVar::runADJ done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFct4DEnsVar<MODEL, OBS>::addIncr(CtrlVar_ & xx, const CtrlInc_ & dx,
                                          PostProcessor<Increment_> &) const {
  Log::trace() << "CostFct4DEnsVar::addIncr start" << std::endl;
  ASSERT(xx.state().validTime() == subWinTime_);
  ASSERT(dx.state().validTime() == subWinTime_);
  xx.state() += dx.state();
  Log::trace() << "CostFct4DEnsVar::addIncr done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFCT4DENSVAR_H_
