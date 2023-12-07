/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020-2021 UCAR.
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
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJbJq.h"
#include "oops/assimilation/CostJcDFI.h"
#include "oops/assimilation/CostJo.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/LinearModel.h"
#include "oops/base/Model.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/State.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
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
  typedef LinearModel<MODEL>              LinearModel_;

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

 protected:
  const Geometry_ & geometry() const override {return *resol_;}

 private:
  void addIncr(CtrlVar_ &, const CtrlInc_ &, PostProcessor<Increment_> &) const override;

  CostJbJq<MODEL, OBS> * newJb(const eckit::Configuration &, const Geometry_ &) const override;
  CostJo<MODEL, OBS>       * newJo(const eckit::Configuration &) const override;
  CostTermBase<MODEL, OBS> * newJc(const eckit::Configuration &, const Geometry_ &) const override;
  void doLinearize(const Geometry_ &, const eckit::Configuration &, CtrlVar_ &, CtrlVar_ &,
                   PostProcessor<State_> &, PostProcessorTLAD<MODEL> &) override;

  size_t nsubwin_;
  size_t nsublocal_;
  const util::TimeWindow timeWindow_;
  util::Duration subWinLength_;
  std::vector<util::DateTime> subWinBgn_;
  std::vector<util::DateTime> subWinEnd_;
  eckit::mpi::Comm * commSpace_;
  eckit::mpi::Comm * commTime_;
  std::unique_ptr<const Geometry_> resol_;
  std::unique_ptr<Model_> model_;
  const Variables ctlvars_;
  std::shared_ptr<LinearModel_> tlm_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostFctWeak<MODEL, OBS>::CostFctWeak(const eckit::Configuration & conf,
                                     const eckit::mpi::Comm & comm)
  : CostFunction<MODEL, OBS>::CostFunction(),
    timeWindow_(conf.getSubConfiguration("time window")),
    subWinLength_(conf.getString("subwindow")),
    subWinBgn_(), subWinEnd_(),
    commSpace_(nullptr), commTime_(nullptr),
    resol_(), model_(), ctlvars_(conf, "analysis variables"), tlm_()
{
  Log::trace() << "CostFctWeak::CostFctWeak start" << std::endl;

  nsubwin_ = timeWindow_.length().toSeconds() / subWinLength_.toSeconds();  // Not like 4D-En-Var
  ASSERT(timeWindow_.length().toSeconds() == subWinLength_.toSeconds()*(int64_t)nsubwin_);

// Define sub-windows
  size_t ntasks = comm.size();
  ASSERT(ntasks % nsubwin_ == 0);
  size_t mysubwin = 0;
  const bool parallel = conf.getBool("parallel subwindows", true);
  if (parallel) {
    nsublocal_ = 1;
    mysubwin = comm.rank() / (ntasks / nsubwin_);
    ASSERT(mysubwin < nsubwin_);
    const util::DateTime mytime = timeWindow_.start() + mysubwin * subWinLength_;
    subWinBgn_.push_back(mytime);
    subWinEnd_.push_back(mytime + subWinLength_);
  } else {
    nsublocal_ = nsubwin_;
    for (size_t jsub = 0; jsub < nsubwin_; ++jsub) {
      const util::DateTime mytime = timeWindow_.start() + jsub * subWinLength_;
      subWinBgn_.push_back(mytime);
      subWinEnd_.push_back(mytime + subWinLength_);
    }
    ASSERT(subWinBgn_[0] == timeWindow_.start());
    ASSERT(subWinEnd_[nsubwin_ - 1] == timeWindow_.end());
    throw eckit::NotImplemented("CostFctWeak: non parallel 4D-En-Var not tested yet", Here());
  }
  ASSERT(subWinBgn_.size() == nsublocal_);
  ASSERT(subWinEnd_.size() == nsublocal_);

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
  ASSERT(commTime_->size() * commSpace_->size() == comm.size());

// Now can setup the rest
  resol_.reset(new Geometry_(eckit::LocalConfiguration(conf, "geometry"),
                             *commSpace_, *commTime_));
  model_.reset(new Model_(*resol_, eckit::LocalConfiguration(conf, "model")));

  this->setupTerms(conf);

  Log::trace() << "CostFctWeak::CostFctWeak done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJbJq<MODEL, OBS> * CostFctWeak<MODEL, OBS>::newJb(const eckit::Configuration & jbConf,
                                                 const Geometry_ & resol) const {
  Log::trace() << "CostFctWeak::newJb start" << std::endl;
  std::vector<util::DateTime> times;
  for (util::DateTime jj = timeWindow_.start(); jj < timeWindow_.end(); jj += subWinLength_) {
    times.push_back(jj);
  }
  return new CostJbJq<MODEL, OBS>(times, jbConf, *commTime_, resol, ctlvars_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJo<MODEL, OBS> * CostFctWeak<MODEL, OBS>::newJo(const eckit::Configuration & joConf) const {
  Log::trace() << "CostFctWeak::newJo" << std::endl;
  return new CostJo<MODEL, OBS>(joConf, *commSpace_,
                                timeWindow_.createSubWindow(subWinBgn_[0],
                                                             subWinEnd_[nsublocal_ - 1]),
                                *commTime_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostTermBase<MODEL, OBS> * CostFctWeak<MODEL, OBS>::newJc(const eckit::Configuration & jcConf,
                                                          const Geometry_ & resol) const {
  Log::trace() << "CostFctWeak::newJc" << std::endl;
  if (nsublocal_ > 1) {
    throw eckit::NotImplemented("CostFctWeak::newJc: no Jc for multiple sub windows", Here());
//  because cannot return more than 1 Jc for now
  }
  const eckit::LocalConfiguration jcdfi(jcConf, "jcdfi");
  const util::DateTime vt(subWinBgn_[0] + subWinLength_/2);
  return new CostJcDFI<MODEL, OBS>(jcdfi, resol, vt, subWinLength_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::runNL(CtrlVar_ & xx, PostProcessor<State_> & post) const {
  Log::trace() << "CostFctWeak::runNL start" << std::endl;
  ASSERT(xx.states().is_4d());
  for (size_t jsub = 0; jsub < nsublocal_; ++jsub) {
    ASSERT(xx.state(jsub).validTime() == subWinBgn_[jsub]);

    model_->forecast(xx.state(jsub), xx.modVar(), subWinLength_, post);

    ASSERT(xx.state().validTime() == subWinEnd_[jsub]);
  }
  Log::trace() << "CostFctWeak::runNL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::doLinearize(const Geometry_ & resol,
                                          const eckit::Configuration & innerConf,
                                          CtrlVar_ & bg, CtrlVar_ & fg,
                                          PostProcessor<State_> & pp,
                                          PostProcessorTLAD<MODEL> & pptraj) {
  Log::trace() << "CostFctWeak::doLinearize start" << std::endl;
  ASSERT(bg.states().is_4d());
  ASSERT(fg.states().is_4d());
// Setup linear model (including trajectory)
  eckit::LocalConfiguration lmConf(innerConf, "linear model");
  tlm_.reset(new LinearModel_(resol, lmConf));
  pp.enrollProcessor(new TrajectorySaver<MODEL>(lmConf, resol, fg.modVar(), tlm_, pptraj));

  Log::trace() << "CostFctWeak::doLinearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::runTLM(CtrlInc_ & dx,
                                     PostProcessorTLAD<MODEL> & cost,
                                     PostProcessor<Increment_> post,
                                     const bool idModel) const {
  Log::trace() << "CostFctWeak: runTLM start" << std::endl;
  ASSERT(dx.states().is_4d());
  for (size_t jsub = 0; jsub < nsublocal_; ++jsub) {
    ASSERT(dx.state(jsub).validTime() == subWinBgn_[jsub]);

    Variables incvars = dx.state(jsub).variables();
    tlm_->forecastTL(dx.state(jsub), dx.modVar(), subWinLength_, post, cost, idModel);

    Log::info() << "CostFctWeak::runTLM: " << dx.state(jsub) << std::endl;
    ASSERT(dx.state(jsub).validTime() == subWinEnd_[jsub]);
  }
  Log::trace() << "CostFctWeak: runTLM done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::runTLM(CtrlInc_ & dx, const bool idModel) const {
  Log::trace() << "CostFctWeak: runTLM start" << std::endl;
  ASSERT(dx.states().is_4d());
  PostProcessor<Increment_> post;
  PostProcessorTLAD<MODEL> cost;

  for (size_t jsub = 0; jsub < nsublocal_; ++jsub) {
    ASSERT(dx.state(jsub).validTime() == subWinBgn_[jsub]);

    if (idModel) {
      dx.state(jsub).updateTime(subWinLength_);
    } else {
      Variables incvars = dx.state(jsub).variables();
      tlm_->forecastTL(dx.state(jsub), dx.modVar(), subWinLength_, post, cost);
    }

    Log::info() << "CostFctWeak::runTLM: " << dx.state(jsub) << std::endl;
    ASSERT(dx.state(jsub).validTime() == subWinEnd_[jsub]);
  }
  Log::trace() << "CostFctWeak: runTLM done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::zeroAD(CtrlInc_ & dx) const {
  Log::trace() << "CostFctWeak::zeroAD start" << std::endl;
  ASSERT(dx.states().is_4d());
  for (size_t jsub = 0; jsub < nsublocal_; ++jsub) {
    dx.state(jsub).zero(subWinEnd_[jsub]);
  }
  dx.modVar().zero();
  dx.obsVar().zero();
  Log::trace() << "CostFctWeak::zeroAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::runADJ(CtrlInc_ & dx,
                                     PostProcessorTLAD<MODEL> & cost,
                                     PostProcessor<Increment_> post,
                                     const bool idModel) const {
  Log::trace() << "CostFctWeak: runADJ start" << std::endl;
  ASSERT(dx.states().is_4d());
  for (size_t jsub = 0; jsub < nsublocal_; ++jsub) {
    ASSERT(dx.state(jsub).validTime() == subWinEnd_[jsub]);

    Variables incvars = dx.state(jsub).variables();
    tlm_->forecastAD(dx.state(jsub), dx.modVar(), subWinLength_, post, cost, idModel);

    Log::info() << "CostFctWeak::runADJ: " << dx.state(jsub) << std::endl;
    ASSERT(dx.state().validTime() == subWinBgn_[jsub]);
  }
  Log::trace() << "CostFctWeak: runADJ done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::runADJ(CtrlInc_ & dx, const bool idModel) const {
  Log::trace() << "CostFctWeak::runADJ start" << std::endl;
  ASSERT(dx.states().is_4d());
  PostProcessor<Increment_> post;
  PostProcessorTLAD<MODEL> cost;

  for (size_t jsub = 0; jsub < nsublocal_; ++jsub) {
    ASSERT(dx.state(jsub).validTime() == subWinEnd_[jsub]);

    if (idModel) {
      dx.state(jsub).updateTime(-subWinLength_);
    } else {
      Variables incvars = dx.state(jsub).variables();
      tlm_->forecastAD(dx.state(jsub), dx.modVar(), subWinLength_, post, cost);
    }

    Log::info() << "CostFctWeak::runADJ: " << dx.state(jsub) << std::endl;
    ASSERT(dx.state().validTime() == subWinBgn_[jsub]);
  }
  Log::trace() << "CostFctWeak::runADJ done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::addIncr(CtrlVar_ & xx, const CtrlInc_ & dx,
                                      PostProcessor<Increment_> & post) const {
  Log::trace() << "CostFctWeak::addIncr start" << std::endl;
  ASSERT(xx.states().is_4d());
  ASSERT(dx.states().is_4d());
  for (size_t jsub = 0; jsub < nsublocal_; ++jsub) {
    ASSERT(xx.state(jsub).validTime() == subWinBgn_[jsub]);
    ASSERT(dx.state(jsub).validTime() == subWinBgn_[jsub]);
  }
  xx.states() += dx.states();
  Log::trace() << "CostFctWeak::addIncr done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFCTWEAK_H_
