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

#ifndef OOPS_ASSIMILATION_COSTFCT4DVAR_H_
#define OOPS_ASSIMILATION_COSTFCT4DVAR_H_

#include <memory>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/CostJb3D.h"
#include "oops/assimilation/CostJcDFI.h"
#include "oops/assimilation/CostJo.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/LinearModel.h"
#include "oops/base/Model.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/State.h"
#include "oops/base/TrajectorySaver.h"
#include "oops/base/VariableChangeParametersBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/LinearVariableChange.h"
#include "oops/interface/VariableChange.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

template <typename MODEL, typename OBS> class CostTermBase;

/// Parameters for the 4D-Var cost function
template <typename MODEL>
class CostFct4DVarParameters : public CostFunctionParametersBase<MODEL> {
  OOPS_CONCRETE_PARAMETERS(CostFct4DVarParameters, CostFunctionParametersBase<MODEL>)

 public:
  typedef typename State<MODEL>::Parameters_           StateParameters_;
  typedef ModelParametersWrapper<MODEL>                ModelParameters_;
  typedef ModelSpaceCovarianceParametersWrapper<MODEL> CovarianceParameters_;

  RequiredParameter<ModelParameters_> model{"model", "model", this};

  // Variable Change
  OptionalParameter<eckit::LocalConfiguration> variableChange{"variable change", this};

  // options for Jb term
  RequiredParameter<StateParameters_> background{"background", "background state", this};
  RequiredParameter<CovarianceParameters_> backgroundError{"background error", "background error",
      this};
  OptionalParameter<eckit::LocalConfiguration> modelAuxControl{"model aux control", this};
  OptionalParameter<eckit::LocalConfiguration> modelAuxError{"model aux error", this};
};

/// Strong Constraint 4D-Var Cost Function
/*!
 * This class is not really necessary since it is only a special
 * case of the more general weak constraint 4D-Var cost function
 * with one sub-window. It is provided for readability.
 */

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS> class CostFct4DVar : public CostFunction<MODEL, OBS> {
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
  typedef CostFct4DVarParameters<MODEL>   Parameters_;

  CostFct4DVar(const Parameters_ &, const eckit::mpi::Comm &);
  ~CostFct4DVar() {}

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

  CostJb3D<MODEL>     * newJb(const eckit::Configuration &, const Geometry_ &,
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
  const eckit::mpi::Comm & comm_;
  Geometry_ resol_;
  Model_ model_;
  const Variables ctlvars_;
  std::shared_ptr<LinearModel_> tlm_;
  VarCha_ an2model_;
  std::unique_ptr<LinVarCha_> inc2model_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostFct4DVar<MODEL, OBS>::CostFct4DVar(const Parameters_ & params,
                                       const eckit::mpi::Comm & comm)
  : CostFunction<MODEL, OBS>::CostFunction(), comm_(comm),
    resol_(params.geometry, comm),
    model_(resol_, params.model.value().modelParameters),
    ctlvars_(params.analysisVariables), tlm_(),
    an2model_(params.variableChange.value() != boost::none ?
              *params.variableChange.value() :
              eckit::LocalConfiguration(), resol_),
    inc2model_()
{
  Log::trace() << "CostFct4DVar:CostFct4DVar" << std::endl;
  windowLength_ = params.windowLength;
  windowBegin_ = params.windowBegin;
  windowEnd_ = windowBegin_ + windowLength_;
  this->setupTerms(params.toConfiguration());
  // ASSERT(ctlvars_ <= this->background().state().variables());
  Log::trace() << "CostFct4DVar constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJb3D<MODEL> * CostFct4DVar<MODEL, OBS>::newJb(const eckit::Configuration & jbConf,
                                                  const Geometry_ & resol,
                                                  const CtrlVar_ & xb) const {
  return new CostJb3D<MODEL>(jbConf, resol, ctlvars_, windowLength_, xb.state());
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostJo<MODEL, OBS> * CostFct4DVar<MODEL, OBS>::newJo(const eckit::Configuration & joConf) const {
  return new CostJo<MODEL, OBS>(joConf, comm_, windowBegin_, windowEnd_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
CostTermBase<MODEL, OBS> * CostFct4DVar<MODEL, OBS>::newJc(const eckit::Configuration & jcConf,
                                                           const Geometry_ & resol) const {
  const eckit::LocalConfiguration jcdfi(jcConf, "jcdfi");
  const util::DateTime vt(windowBegin_ + windowLength_/2);
  return new CostJcDFI<MODEL, OBS>(jcdfi, resol, vt, windowLength_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct4DVar<MODEL, OBS>::runNL(CtrlVar_ & xx, PostProcessor<State_> & post) const {
  ASSERT(xx.state().validTime() == windowBegin_);

  Variables anvars(xx.state().variables());
  an2model_.changeVar(xx.state(), model_.variables());
  model_.forecast(xx.state(), xx.modVar(), windowLength_, post);
  an2model_.changeVarInverse(xx.state(), anvars);
  ASSERT(xx.state().validTime() == windowEnd_);
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFct4DVar<MODEL, OBS>::doLinearize(const Geometry_ & resol,
                                           const eckit::Configuration & innerConf,
                                           const CtrlVar_ & bg, const CtrlVar_ & fg,
                                           PostProcessor<State_> & pp,
                                           PostProcessorTLAD<MODEL> & pptraj) {
  Log::trace() << "CostFct4DVar::doLinearize start" << std::endl;
  eckit::LocalConfiguration lmConf(innerConf, "linear model");
  // Setup linear model (and trajectory)
  tlm_.reset(new LinearModel_(resol, lmConf));
  pp.enrollProcessor(new TrajectorySaver<MODEL>(lmConf, resol, fg.modVar(), tlm_, pptraj));

  // Create variable change
  inc2model_.reset(new LinVarCha_(resol, innerConf.getSubConfiguration("linear variable change")));

  // Trajectory for linear variable change
  inc2model_->setTrajectory(bg.state(), fg.state());

  Log::trace() << "CostFct4DVar::doLinearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct4DVar<MODEL, OBS>::runTLM(CtrlInc_ & dx,
                                      PostProcessorTLAD<MODEL> & cost,
                                      PostProcessor<Increment_> post,
                                      const bool idModel) const {
  Log::trace() << "CostFct4DVar::runTLM starting" << std::endl;
  ASSERT(dx.state().validTime() == windowBegin_);

  Variables incvars = dx.state().variables();

  inc2model_->multiply(dx.state(), tlm_->variables());
  tlm_->forecastTL(dx.state(), dx.modVar(), windowLength_, post, cost, idModel);
  inc2model_->multiplyInverse(dx.state(), incvars);
  ASSERT(dx.state().validTime() == windowEnd_);
  Log::trace() << "CostFct4DVar::runTLM done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct4DVar<MODEL, OBS>::zeroAD(CtrlInc_ & dx) const {
  dx.state().zero(windowEnd_);
  dx.modVar().zero();
  dx.obsVar().zero();
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void CostFct4DVar<MODEL, OBS>::runADJ(CtrlInc_ & dx,
                                      PostProcessorTLAD<MODEL> & cost,
                                      PostProcessor<Increment_> post,
                                      const bool idModel) const {
  Log::trace() << "CostFct4DVar::runADJ starting" << std::endl;
  ASSERT(dx.state().validTime() == windowEnd_);

  Variables incvars = dx.state().variables();
  inc2model_->multiplyInverseAD(dx.state(), tlm_->variables());
  tlm_->forecastAD(dx.state(), dx.modVar(), windowLength_, post, cost, idModel);
  inc2model_->multiplyAD(dx.state(), incvars);

  ASSERT(dx.state().validTime() == windowBegin_);
  Log::trace() << "CostFct4DVar::runADJ done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostFct4DVar<MODEL, OBS>::addIncr(CtrlVar_ & xx, const CtrlInc_ & dx,
                                       PostProcessor<Increment_> &) const {
  xx.state() += dx.state();
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTFCT4DVAR_H_
