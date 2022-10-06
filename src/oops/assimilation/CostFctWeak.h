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
#include "oops/base/StateParametersND.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/base/Variables.h"
#include "oops/interface/LinearVariableChange.h"
#include "oops/interface/VariableChange.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

template <typename MODEL, typename OBS> class CostTermBase;

/// Parameters for the Weak cost function
template <typename MODEL, typename OBS>
class CostFctWeakParameters : public CostFunctionParametersBase<MODEL, OBS> {
  // This typedef prevents the macro below from choking on the 2 args of the templated type
  typedef CostFunctionParametersBase<MODEL, OBS> CostFuntionParametersBase_;
  OOPS_CONCRETE_PARAMETERS(CostFctWeakParameters, CostFuntionParametersBase_);

 public:
  typedef ModelParametersWrapper<MODEL> ModelParameters_;
  typedef StateParameters4D<MODEL>      StateParameters4D_;
  typedef typename VariableChange<MODEL>::Parameters_  VariableChangeParameters_;

  RequiredParameter<ModelParameters_> model{"model", "model", this};
  RequiredParameter<util::Duration> subwindow{"subwindow", "length of assimilation subwindows",
      this};

  Parameter<VariableChangeParameters_> variableChange{"variable change",
           "variable change from B matrix variables to model variables", {}, this};

  // options for Jb term
  RequiredParameter<StateParameters4D_> background{"background", "background state(s)", this};
  // Currently `ModelSpaceCovarianceParametersWrapper` doesn't support multiple covariances for
  // multiple models, so read this as a config.
  RequiredParameter<eckit::LocalConfiguration> backgroundError{"background error",
      "background error(s)", this};
};

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
  typedef VariableChange<MODEL>           VarCha_;
  typedef LinearVariableChange<MODEL>     LinVarCha_;

 public:
  typedef CostFctWeakParameters<MODEL, OBS> Parameters_;

  CostFctWeak(const Parameters_ &, const eckit::mpi::Comm &);
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
  CostJo<MODEL, OBS>       * newJo(const ObserversParameters<MODEL, OBS> &) const override;
  CostTermBase<MODEL, OBS> * newJc(const eckit::Configuration &, const Geometry_ &) const override;
  void doLinearize(const Geometry_ &, const eckit::Configuration &,
                   const CtrlVar_ &, const CtrlVar_ &,
                   PostProcessor<State_> &, PostProcessorTLAD<MODEL> &) override;
  const Geometry_ & geometry() const override {return *resol_;}

  util::Duration subWinLength_;
  util::DateTime subWinBegin_;
  util::DateTime subWinEnd_;
  size_t nsubwin_;
  size_t mysubwin_;
  std::unique_ptr<const Geometry_> resol_;
  std::unique_ptr<Model_> model_;
  const Variables ctlvars_;
  std::shared_ptr<LinearModel_> tlm_;
  std::unique_ptr<VarCha_> an2model_;
  std::unique_ptr<LinVarCha_> inc2model_;
  eckit::mpi::Comm * commSpace_;
  eckit::mpi::Comm * commTime_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostFctWeak<MODEL, OBS>::CostFctWeak(const Parameters_ & params,
                                     const eckit::mpi::Comm & comm)
  : CostFunction<MODEL, OBS>::CostFunction(), resol_(), model_(),
    ctlvars_(params.analysisVariables), tlm_(), an2model_(),
    inc2model_(), commSpace_(nullptr), commTime_(nullptr)
{
  const util::Duration windowLength = params.windowLength;
  const util::DateTime windowBegin = params.windowBegin;
  subWinLength_ = params.subwindow;

  nsubwin_ = windowLength.toSeconds() / subWinLength_.toSeconds();
  ASSERT(windowLength.toSeconds() == subWinLength_.toSeconds()*(int64_t)nsubwin_);

  size_t ntasks = comm.size();
  ASSERT(ntasks % nsubwin_ == 0);
  size_t myrank = comm.rank();
  size_t ntaskpslot = ntasks / nsubwin_;
  size_t mysubwin_ = myrank / ntaskpslot;

// Define local sub-window
  subWinBegin_ = windowBegin + mysubwin_ * subWinLength_;
  subWinEnd_ = subWinBegin_ + subWinLength_;

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
  resol_.reset(new Geometry_(params.geometry, *commSpace_, *commTime_));
  model_.reset(new Model_(*resol_, params.model.value().modelParameters));
  an2model_.reset(new VarCha_(params.variableChange, *resol_));
  this->setupTerms(params.toConfiguration());

  Log::trace() << "CostFctWeak constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJbJq<MODEL> * CostFctWeak<MODEL, OBS>::newJb(const eckit::Configuration & jbConf,
                                                 const Geometry_ & resol,
                                                 const CtrlVar_ & xb) const {
  return new CostJbJq<MODEL>(jbConf, *commTime_, resol, ctlvars_, xb.state());
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJo<MODEL, OBS> * CostFctWeak<MODEL, OBS>::newJo(
    const ObserversParameters<MODEL, OBS> & joParams) const {
  return new CostJo<MODEL, OBS>(joParams, *commSpace_,
                                subWinBegin_, subWinEnd_, *commTime_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostTermBase<MODEL, OBS> * CostFctWeak<MODEL, OBS>::newJc(const eckit::Configuration & jcConf,
                                                          const Geometry_ & resol) const {
  const eckit::LocalConfiguration jcdfi(jcConf, "jcdfi");
  const util::DateTime vt(subWinBegin_ + subWinLength_/2);
  return new CostJcDFI<MODEL, OBS>(jcdfi, resol, vt, subWinLength_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::runNL(CtrlVar_ & xx, PostProcessor<State_> & post) const {
  ASSERT(xx.state().validTime() == subWinBegin_);

  Variables anvars(xx.state().variables());
  an2model_->changeVar(xx.state(), model_->variables());
  model_->forecast(xx.state(), xx.modVar(), subWinLength_, post);
  an2model_->changeVarInverse(xx.state(), anvars);

  ASSERT(xx.state().validTime() == subWinEnd_);
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::doLinearize(const Geometry_ & resol,
                                          const eckit::Configuration & innerConf,
                                          const CtrlVar_ & bg, const CtrlVar_ & fg,
                                          PostProcessor<State_> & pp,
                                          PostProcessorTLAD<MODEL> & pptraj) {
  Log::trace() << "CostFctWeak::doLinearize start" << std::endl;
  eckit::LocalConfiguration lmConf(innerConf, "linear model");
  // Setup linear model (and trajectory)
  tlm_.reset(new LinearModel_(resol, lmConf));
  pp.enrollProcessor(new TrajectorySaver<MODEL>(lmConf, resol, fg.modVar(), tlm_, pptraj));

  // Create variable change
  std::unique_ptr<eckit::LocalConfiguration> lvcConf;
  inc2model_.reset(new LinVarCha_(resol, innerConf.getSubConfiguration("linear variable change")));

  // Trajecotry for linear variable change
  inc2model_->changeVarTraj(fg.state(), tlm_->variables());

  Log::trace() << "CostFctWeak::doLinearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::runTLM(CtrlInc_ & dx,
                                     PostProcessorTLAD<MODEL> & cost,
                                     PostProcessor<Increment_> post,
                                     const bool idModel) const {
  Log::trace() << "CostFctWeak: runTLM start" << std::endl;

  ASSERT(dx.state().validTime() == subWinBegin_);

  Variables incvars = dx.state().variables();

  inc2model_->changeVarTL(dx.state(), tlm_->variables());
  tlm_->forecastTL(dx.state(), dx.modVar(), subWinLength_, post, cost, idModel);
  inc2model_->changeVarInverseTL(dx.state(), incvars);

  ASSERT(dx.state().validTime() == subWinEnd_);
  Log::trace() << "CostFctWeak: runTLM done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::runTLM(CtrlInc_ & dx, const bool idModel) const {
  PostProcessor<Increment_> post;
  PostProcessorTLAD<MODEL> cost;

  ASSERT(dx.state().validTime() == subWinBegin_);

  Variables incvars = dx.state().variables();
  if (idModel) {
    dx.state().updateTime(subWinLength_);
  } else {
    inc2model_->changeVarTL(dx.state(), tlm_->variables());
    tlm_->forecastTL(dx.state(), dx.modVar(), subWinLength_, post, cost);
    inc2model_->changeVarInverseTL(dx.state(), incvars);
  }

  ASSERT(dx.state().validTime() == subWinEnd_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::zeroAD(CtrlInc_ & dx) const {
  dx.state().zero(subWinEnd_);
  dx.modVar().zero();
  dx.obsVar().zero();
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::runADJ(CtrlInc_ & dx,
                                     PostProcessorTLAD<MODEL> & cost,
                                     PostProcessor<Increment_> post,
                                     const bool idModel) const {
  Log::trace() << "CostFctWeak: runADJ start" << std::endl;

  ASSERT(dx.state().validTime() == subWinEnd_);

  Variables incvars = dx.state().variables();
  inc2model_->changeVarInverseAD(dx.state(), tlm_->variables());
  tlm_->forecastAD(dx.state(), dx.modVar(), subWinLength_, post, cost, idModel);
  inc2model_->changeVarAD(dx.state(), incvars);

  ASSERT(dx.state().validTime() == subWinBegin_);
  Log::trace() << "CostFctWeak: runADJ done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::runADJ(CtrlInc_ & dx, const bool idModel) const {
  PostProcessor<Increment_> post;
  PostProcessorTLAD<MODEL> cost;

  ASSERT(dx.state().validTime() == subWinEnd_);

  Variables incvars = dx.state().variables();
  if (idModel) {
    dx.state().updateTime(-subWinLength_);
  } else {
    inc2model_->changeVarInverseAD(dx.state(), tlm_->variables());
    tlm_->forecastAD(dx.state(), dx.modVar(), subWinLength_, post, cost);
    inc2model_->changeVarAD(dx.state(), incvars);
  }

  ASSERT(dx.state().validTime() == subWinBegin_);
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFctWeak<MODEL, OBS>::addIncr(CtrlVar_ & xx, const CtrlInc_ & dx,
                                      PostProcessor<Increment_> & post) const {
  xx.state() += dx.state();
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFCTWEAK_H_
