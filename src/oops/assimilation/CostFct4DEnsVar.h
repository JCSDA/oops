/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2023 UCAR
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
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJb4D.h"
#include "oops/assimilation/CostJcDFI.h"
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

  void runTLM(CtrlInc_ &, PostProcessorTLAD<MODEL> &, PostProcessor<Increment_>,
              const bool idModel = false) const override;
  void runADJ(CtrlInc_ &, PostProcessorTLAD<MODEL> &, PostProcessor<Increment_>,
              const bool idModel = false) const override;
  void zeroAD(CtrlInc_ &) const override;

  void runNL(CtrlVar_ &, PostProcessor<State_>&) const override;

 protected:
  const Geometry_ & geometry() const override {return *resol_;}

 private:
  void addIncr(CtrlVar_ &, const CtrlInc_ &, PostProcessor<Increment_>&) const override;

  CostJb4D<MODEL, OBS> * newJb(const eckit::Configuration &, const Geometry_ &) const override;
  CostJo<MODEL, OBS>       * newJo(const eckit::Configuration &) const override;
  CostTermBase<MODEL, OBS> * newJc(const eckit::Configuration &, const Geometry_ &) const override;
  void doLinearize(const Geometry_ &, const eckit::Configuration &, CtrlVar_ &, CtrlVar_ &,
                   PostProcessor<State_> &, PostProcessorTLAD<MODEL> &) override;

  size_t nsubwin_;
  size_t nsublocal_;
  size_t last_;
  const util::TimeWindow timeWindow_;
  util::Duration subWinLength_;
  std::vector<util::DateTime> subWinTime_;
  std::vector<util::DateTime> subWinBgn_;
  std::vector<util::DateTime> subWinEnd_;
  eckit::mpi::Comm * commSpace_;
  eckit::mpi::Comm * commTime_;
  std::unique_ptr<const Geometry_> resol_;
  const Variables ctlvars_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostFct4DEnsVar<MODEL, OBS>::CostFct4DEnsVar(const eckit::Configuration & conf,
                                             const eckit::mpi::Comm & comm)
  : CostFunction<MODEL, OBS>::CostFunction(),
    timeWindow_(conf.getSubConfiguration("time window")),
    subWinLength_(conf.getString("subwindow")),
    subWinTime_(), subWinBgn_(), subWinEnd_(),
    resol_(), ctlvars_(conf, "analysis variables")
{
  Log::trace() << "CostFct4DEnsVar::CostFct4DEnsVar start" << std::endl;

  nsubwin_ = timeWindow_.length().toSeconds() / subWinLength_.toSeconds() + 1;  // Not like WC
  ASSERT(timeWindow_.length().toSeconds() == subWinLength_.toSeconds() * (int64_t)(nsubwin_ - 1));

// Define sub-windows
  size_t ntasks = comm.size();
  size_t mysubwin = 0;
  const bool parallel = conf.getBool("parallel subwindows", true);
  if (parallel) {
    ASSERT(ntasks % nsubwin_ == 0);
    nsublocal_ = 1;
    mysubwin = comm.rank() / (ntasks / nsubwin_);
    ASSERT(mysubwin < nsubwin_);
    const util::DateTime mytime = timeWindow_.start() + mysubwin * subWinLength_;
    subWinTime_.push_back(mytime);
    subWinBgn_.push_back(mytime - subWinLength_/2);
    subWinEnd_.push_back(mytime + subWinLength_/2);
    if (mysubwin == 0) subWinBgn_[0] = timeWindow_.start();
    if (mysubwin == nsubwin_ - 1) subWinEnd_[0] = timeWindow_.end();
  } else {
    nsublocal_ = nsubwin_;
    for (size_t jsub = 0; jsub < nsubwin_; ++jsub) {
      const util::DateTime mytime = timeWindow_.start() + jsub * subWinLength_;
      subWinTime_.push_back(mytime);
      subWinBgn_.push_back(mytime - subWinLength_/2);
      subWinEnd_.push_back(mytime + subWinLength_/2);
    }
    subWinBgn_[0] = timeWindow_.start();
    subWinEnd_[nsubwin_ - 1] = timeWindow_.end();
  }
  last_ = nsublocal_ - 1;

// Create a communicator for same sub-window, to be used for communications in space
  std::string sgeom = "comm_geom_" + std::to_string(mysubwin);
  char const *geomName = sgeom.c_str();
  commSpace_ = &comm.split(mysubwin, geomName);

// Create a communicator for same local area, to be used for communications in time
  size_t myarea = commSpace_->rank();
  std::string stime = "comm_time_" + std::to_string(myarea);
  char const *timeName = stime.c_str();
  commTime_ = &comm.split(myarea, timeName);
  ASSERT(commTime_->size() == (nsubwin_ / nsublocal_));

// Now can setup the rest
  resol_.reset(new Geometry_(eckit::LocalConfiguration(conf, "geometry"),
                             *commSpace_, *commTime_));

  this->setupTerms(conf);

  Log::trace() << "CostFct4DEnsVar::CostFct4DEnsVar done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJb4D<MODEL, OBS> * CostFct4DEnsVar<MODEL, OBS>::newJb(const eckit::Configuration & jbConf,
                                                          const Geometry_ & resol) const {
  Log::trace() << "CostFct4DEnsVar::newJb" << std::endl;
  std::vector<util::DateTime> times;
  for (util::DateTime jj = timeWindow_.start(); jj <= timeWindow_.end(); jj += subWinLength_) {
    times.push_back(jj);
  }
  return new CostJb4D<MODEL, OBS>(times, jbConf, *commTime_, resol, ctlvars_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJo<MODEL, OBS> * CostFct4DEnsVar<MODEL, OBS>::newJo(const eckit::Configuration & joConf) const {
  Log::trace() << "CostFct4DEnsVar::newJo" << std::endl;
  return new CostJo<MODEL, OBS>(joConf, *commSpace_,
                                timeWindow_.createSubWindow(subWinBgn_[0], subWinEnd_[last_]),
                                *commTime_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostTermBase<MODEL, OBS> * CostFct4DEnsVar<MODEL, OBS>::newJc(const eckit::Configuration & jcConf,
                                                              const Geometry_ & resol) const {
  Log::trace() << "CostFct4DEnsVar::newJc" << std::endl;
  if (nsublocal_ != nsubwin_) {
    throw eckit::NotImplemented("CostFct4DEnsVar::newJc: no parallel Jc", Here());
  }
  const eckit::LocalConfiguration jcdfi(jcConf, "jcdfi");
  const util::DateTime vt(timeWindow_.start() + timeWindow_.length()/2);
  return new CostJcDFI<MODEL, OBS>(jcdfi, resol, vt, timeWindow_.length(), subWinLength_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct4DEnsVar<MODEL, OBS>::runNL(CtrlVar_ & xx, PostProcessor<State_> & post) const {
  Log::trace() << "CostFct4DEnsVar::runNL start" << std::endl;
  ASSERT(xx.states().is_4d());
  post.initialize(xx.state(0), subWinEnd_[last_], subWinLength_);
  for (size_t jsub = 0; jsub < nsublocal_; ++jsub) {
    ASSERT(xx.state(jsub).validTime() == subWinTime_[jsub]);
    post.process(xx.state(jsub));
    Log::info() << "CostFct4DEnsVar::runNL: " << xx.state(jsub)  << std::endl;
  }
  post.finalize(xx.state(last_));
  Log::trace() << "CostFct4DEnsVar::runNL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFct4DEnsVar<MODEL, OBS>::doLinearize(const Geometry_ & resol,
                                              const eckit::Configuration & conf,
                                              CtrlVar_ &, CtrlVar_ &,
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
  ASSERT(dx.states().is_4d());

  cost.initializeTL(dx.state(0), subWinTime_[last_], subWinLength_);
  post.initialize(dx.state(0), subWinTime_[last_], subWinLength_);

  for (size_t jsub = 0; jsub < nsublocal_; ++jsub) {
    ASSERT(dx.state(jsub).validTime() == subWinTime_[jsub]);

    cost.processTL(dx.state(jsub));
    post.process(dx.state(jsub));

    Log::info() << "CostFct4DEnsVar::runTLM: " << dx.state(jsub) << std::endl;
  }

  cost.finalizeTL(dx.state(last_));
  post.finalize(dx.state(last_));

  Log::trace() << "CostFct4DEnsVar::runTLM done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct4DEnsVar<MODEL, OBS>::zeroAD(CtrlInc_ & dx) const {
  Log::trace() << "CostFct4DEnsVar::zeroAD start" << std::endl;
  ASSERT(dx.states().is_4d());
  for (size_t jsub = 0; jsub < nsublocal_; ++jsub) {
    dx.state(jsub).zero(subWinTime_[jsub]);
  }
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
  ASSERT(dx.states().is_4d());
  ASSERT(dx.states().local_time_size() == nsublocal_);

  post.initialize(dx.state(0), subWinTime_[last_], subWinLength_);
  cost.initializeAD(dx.state(0), subWinTime_[last_], subWinLength_);

  for (size_t jsub = 0; jsub < nsublocal_; ++jsub) {
    ASSERT(dx.state(jsub).validTime() == subWinTime_[jsub]);

    cost.processAD(dx.state(jsub));
    post.process(dx.state(jsub));

    Log::info() << "CostFct4DEnsVar::runADJ: " << dx.state(jsub) << std::endl;
  }

  cost.finalizeAD(dx.state(last_));
  post.finalize(dx.state(last_));

  Log::trace() << "CostFct4DEnsVar::runADJ done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFct4DEnsVar<MODEL, OBS>::addIncr(CtrlVar_ & xx, const CtrlInc_ & dx,
                                          PostProcessor<Increment_> &) const {
  Log::trace() << "CostFct4DEnsVar::addIncr start" << std::endl;
  ASSERT(xx.states().is_4d());
  ASSERT(dx.states().is_4d());
  for (size_t jsub = 0; jsub < nsublocal_; ++jsub) {
    ASSERT(xx.state(jsub).validTime() == subWinTime_[jsub]);
    ASSERT(dx.state(jsub).validTime() == subWinTime_[jsub]);
  }
  xx.states() += dx.states();
  Log::trace() << "CostFct4DEnsVar::addIncr done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFCT4DENSVAR_H_
