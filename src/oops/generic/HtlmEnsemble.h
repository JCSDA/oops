/*
 * (C) Copyright 2022-2023 UCAR.
 * (C) Crown copyright 2022-2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_HTLMENSEMBLE_H_
#define OOPS_GENERIC_HTLMENSEMBLE_H_

#include <memory>
#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/Increment4D.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/Model.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/generic/HtlmSimplifiedLinearModel.h"

namespace oops {

template <typename MODEL>
class StatePerturbationParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(StatePerturbationParameters, Parameters)
  typedef ModelSpaceCovarianceParametersWrapper<MODEL> CovarianceParameters_;

 public:
  RequiredParameter<CovarianceParameters_> backgroundError{"background error", this};
  RequiredParameter<Variables> variables{"variables", this};
};

template<typename MODEL>
class NonLinearEnsembleParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(NonLinearEnsembleParameters, Parameters)
  typedef StatePerturbationParameters<MODEL>    StatePerturbationParameters_;
  typedef StateEnsembleParameters<MODEL>        StateEnsembleParameters_;

 public:
    // Configurations of all perturbed members if reading
    OptionalParameter<StateEnsembleParameters_> statesReadIn{"perturbed members",
                                 "read perturbed members from disk", this};
    // Parameters to generate initial ensemble if not reading in
    OptionalParameter<StatePerturbationParameters_>
    statesPerturbation{"ensemble perturbation", "generate ensemble from error covariance", this};

    void check() const;
};

template <typename MODEL>
void NonLinearEnsembleParameters<MODEL>::check() const
{
  if (statesReadIn.value() == boost::none && statesPerturbation.value() == boost::none) {
    ABORT("HtlmEnsembleParameters: both initial ensemble and ensemble perturbation are missing");
  }
  if (statesReadIn.value() != boost::none && statesPerturbation.value() != boost::none) {
    ABORT("HtlmEnsembleParameters: both initial ensemble and ensemble perturbation are present");
  }
}

template <typename MODEL>
class HtlmEnsembleParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(HtlmEnsembleParameters, Parameters);
  typedef typename State<MODEL>::Parameters_    StateParameters_;
  typedef StateEnsembleParameters<MODEL>        StateEnsembleParameters_;
  typedef typename Geometry<MODEL>::Parameters_ GeometryParameters_;
  typedef NonLinearEnsembleParameters<MODEL>    NlEnsParameters_;

 public:
  // Nonlinear forecast model.
  RequiredParameter<eckit::LocalConfiguration> model{"model", this};
  // Background and analysis geometry.
  RequiredParameter<GeometryParameters_> stateGeometry{"state geometry", this};
  // Augmented model state
  Parameter<eckit::LocalConfiguration> modelAuxControl{"model aux control",
    eckit::LocalConfiguration(), this};
  // Augmented Model Increment
  Parameter<eckit::LocalConfiguration> modelAuxIncrement{"model aux increment",
    eckit::LocalConfiguration(), this};
  // Configurations of the control member.
  RequiredParameter<StateParameters_> control{"control member", this};
  // Number of perturbed ensemble members
  RequiredParameter<size_t> ensembleSize{"ensemble size", this};
  // Nonlinear ensemble initalization parameters
  RequiredParameter<NlEnsParameters_> nlEnsemble{"non linear ensemble", this};
};

template <typename MODEL>
class HtlmEnsemble{
  typedef Geometry<MODEL>                               Geometry_;
  typedef Model<MODEL>                                  Model_;
  typedef ModelAuxControl<MODEL>                        ModelAux_;
  typedef ModelAuxIncrement<MODEL>                      ModelAuxIncrement_;
  typedef State<MODEL>                                  State_;
  typedef State4D<MODEL>                                State4D_;
  typedef StateEnsemble<MODEL>                          StateEnsemble_;
  typedef Increment<MODEL>                              Increment_;
  typedef Increment4D<MODEL>                            Increment4D_;
  typedef IncrementEnsemble<MODEL>                      IncrementEnsemble_;
  typedef HtlmEnsembleParameters<MODEL>                 HtlmEnsembleParameters_;
  typedef ModelSpaceCovarianceBase<MODEL>               CovarianceBase_;
  typedef CovarianceFactory<MODEL>                      CovarianceFactory_;
  typedef ModelSpaceCovarianceParametersBase<MODEL>     CovarianceParametersBase_;
  typedef HtlmSimplifiedLinearModel<MODEL>              HtlmSimplifiedLinearModel_;

 public:
  static const std::string classname() {return "oops::HtlmEnsemble";}

  HtlmEnsemble(const HtlmEnsembleParameters_ &, const Variables &, const Geometry_ &);
  void step(const util::Duration &, HtlmSimplifiedLinearModel_ &);

  IncrementEnsemble_ & getLinearEns() {return linearEnsemble_;}
  IncrementEnsemble_ & getLinearErrDe() {return linearErrorDe_;}
  const size_t size() const {return linearEnsemble_.size();}

 private:
  const size_t ensembleSize_;
  const Geometry_ stateGeometry_;
  const Geometry_ & incrementGeometry_;
  //  Model
  const Model_ model_;
  //  control member IC
  State4D_ controlState_;
  //  Augmented state
  ModelAux_ moderr_;
  // Augmented increment
  ModelAuxIncrement_ modauxinc_;
  //  perturbed member ICs
  StateEnsemble_ perturbedStates_;
  //  Linear Ensemble
  IncrementEnsemble_ linearEnsemble_;
  // Nonlinear Differences
  IncrementEnsemble_ nonLinearDifferences_;
  // linear error (dE in the HTLM paper)
  IncrementEnsemble_ linearErrorDe_;
};

//----------------------------------------------------------------------------------

template<typename MODEL>
HtlmEnsemble<MODEL>::HtlmEnsemble(const HtlmEnsembleParameters_ & params,
                                  const Variables & simplifiedLinearModelVars,
                                  const Geometry_ & updateGeometry)
    : ensembleSize_(params.ensembleSize),
      stateGeometry_(params.stateGeometry.value(), updateGeometry.getComm()),
      incrementGeometry_(updateGeometry),
      model_(stateGeometry_, eckit::LocalConfiguration(params.toConfiguration(), "model")),
      controlState_(stateGeometry_, params.control.value().toConfiguration()),
      moderr_(stateGeometry_, params.modelAuxControl.value()),
      modauxinc_(incrementGeometry_, params.modelAuxIncrement.value()),
      perturbedStates_(params.nlEnsemble.value().statesReadIn.value() != boost::none ?
             StateEnsemble_(stateGeometry_, *params.nlEnsemble.value().statesReadIn.value()) :
                                                   StateEnsemble_(controlState_[0], ensembleSize_)),
      linearEnsemble_(incrementGeometry_, simplifiedLinearModelVars,
                      controlState_[0].validTime(), ensembleSize_),
      nonLinearDifferences_(incrementGeometry_, simplifiedLinearModelVars,
                            controlState_[0].validTime(), ensembleSize_),
      linearErrorDe_(nonLinearDifferences_) {
    Log::trace() << "HtlmEnsemble<MODEL>::HtlmEnsemble() starting"
                 << std::endl;
    // Check ensemble initilization params
    params.nlEnsemble.value().check();

  // If chosesn, initialize nonlinear ensemble from background error.
  if (params.nlEnsemble.value().statesPerturbation.value() != boost::none) {
    // set up variables to perturb
    const Variables vars_(params.nlEnsemble.value().statesPerturbation.value()->variables);
    // set up Bmatrix
    const CovarianceParametersBase_ &covarParams =
    params.nlEnsemble.value().statesPerturbation.value()->
                         backgroundError.value().covarianceParameters;

    std::unique_ptr<CovarianceBase_> Bmat(CovarianceFactory_::create(
                         stateGeometry_, vars_, covarParams, controlState_, controlState_));
    //  Generate perturbed states
    Increment4D_ dx(stateGeometry_, vars_, controlState_.times());
    for (size_t jm = 0; jm < ensembleSize_; ++jm) {
      //  Generate linear ensemble increments
      Bmat->randomize(dx);
      //  Add to control state
      perturbedStates_[jm] += dx[0];
    }
  }

  //  Set up linearEnsemble_
  for (size_t m = 0; m < ensembleSize_; ++m) {
    linearEnsemble_[m].diff(controlState_[0], perturbedStates_[m]);
  }

  Log::trace() << "HtlmEnsemble<MODEL>::HtlmEnsemble() done" << std::endl;
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HtlmEnsemble<MODEL>::step(const util::Duration & tstep,
                               HtlmSimplifiedLinearModel_ & simplifiedLinearModel)  {
    Log::trace() << "HtlmEnsemble<MODEL>::step() starting" << std::endl;
    State_ downsampled_Control(incrementGeometry_, controlState_[0]);
    PostProcessor<State_> post;
    for (util::Duration t(0); t < tstep; t+= simplifiedLinearModel.timeResolution()) {
      simplifiedLinearModel.setSimplifiedTrajectory(controlState_[0], downsampled_Control, moderr_);
      model_.forecast(controlState_[0], moderr_, simplifiedLinearModel.timeResolution(), post);
    }
    for (size_t m = 0; m < ensembleSize_; m++) {
      model_.forecast(perturbedStates_[m], moderr_, tstep, post);
      simplifiedLinearModel.forecastSimplifiedTL(linearEnsemble_[m], modauxinc_, tstep);
    }
    for (size_t m = 0; m < ensembleSize_; ++m) {
        nonLinearDifferences_[m].updateTime(tstep);
        nonLinearDifferences_[m].diff(controlState_[0], (perturbedStates_[m]));
    }
    linearErrorDe_ = nonLinearDifferences_;
    for (size_t m = 0; m < ensembleSize_; ++m) {
        linearErrorDe_[m] -= linearEnsemble_[m];
    }
    Log::trace() << "HtlmEnsemble<MODEL>::step() done" << std::endl;
}

//------------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_HTLMENSEMBLE_H_
