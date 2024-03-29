/*
 * (C) Copyright 2022 MetOffice.
 * (C) Copyright 2021-2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_HTLMENSEMBLE_H_
#define OOPS_GENERIC_HTLMENSEMBLE_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/Model.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/generic/LinearModelBase.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Expect.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"


namespace oops {

// parameters for forwarding the ensemble.
template <typename MODEL>
class HtlmEnsembleParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(HtlmEnsembleParameters, Parameters);
  typedef typename State<MODEL>::Parameters_    StateParameters_;
  typedef ModelParametersWrapper<MODEL>         ModelParameters_;
  typedef typename Geometry<MODEL>::Parameters_ GeometryParameters_;

 public:
  // Nonlinear forecast model.
  RequiredParameter<ModelParameters_> model{"model", this};
  // Background and analysis geometry.
  RequiredParameter<GeometryParameters_> stateGeometry{"state geometry", this};
  // Augmented model state
  Parameter<eckit::LocalConfiguration> modelAuxControl{"model aux control",
    eckit::LocalConfiguration(), this};
  // Linear model.
  RequiredParameter<eckit::LocalConfiguration> linearModel{"coeff linear model", this};
  // Augmented Model Increment
  Parameter<eckit::LocalConfiguration> modelAuxIncrement{"model aux increment",
    eckit::LocalConfiguration(), this};
  // Configurations of the control member.
  RequiredParameter<StateParameters_> control{"control member", this};
  // Configurations of all perturbed members.
  RequiredParameter<std::vector<StateParameters_>>
        perturbedMembers{"perturbed members", this};
  // Number of perturbed ensemble members
  RequiredParameter<size_t> ensembleSize{"ensemble size", this};
};

// class declarations for htlmensemble
template <typename MODEL>
class HtlmEnsemble{
  typedef Geometry<MODEL>                               Geometry_;
  typedef Model<MODEL>                                  Model_;
  typedef LinearModel<MODEL>                            LinearModel_;
  typedef ModelAuxControl<MODEL>                        ModelAux_;
  typedef ModelAuxIncrement<MODEL>                      ModelAuxIncrement_;
  typedef State<MODEL>                                  State_;
  typedef Increment<MODEL>                              Increment_;
  typedef HtlmEnsembleParameters<MODEL>                 HtlmEnsembleParameters_;

 public:
  static const std::string classname() {return "oops::HtlmEnsemble";}


// constructor

  /* The geometry for the state is a yaml parameter
   The geometry passed in is used for the TLM geometry
   This is ultimately first constructed in HybridLinearModel from its parameters */

  HtlmEnsemble(const HtlmEnsembleParameters_ &, const Geometry_ &);


  void step(const util::Duration &);


// Needed accessors for passing info to calculator
  std::vector<Increment_> & getLinearEns() {return linearEnsemble_;}
  std::vector<Increment_> & getLinearErrDe() {return linearErrorDe_;}

 private:
  const HtlmEnsembleParameters_ params_;
  const size_t ensembleSize_;
  const Geometry_ stateGeometry_;
  const Geometry_ & incrementGeometry_;

//  Model
  const Model_ model_;
// ptr to linear model
  LinearModel_ simpleLinearModel_;
//  control member IC
  State_ controlState_;
//  perturbed member ICs
  std::vector<State_> perturbedStates_;
//  Linear Ensemble
  std::vector<Increment_> linearEnsemble_;
// Nonlinear Differences
  std::vector<Increment_> nonLinearDifferences_;
// linear error (dE in the HTLM paper)
  std::vector<Increment_> linearErrorDe_;
//  Augmented state
  ModelAux_ moderr_;
// Augmented increment
  ModelAuxIncrement_ modauxinc_;
};

//----------------------------------------------------------------------------------

template<typename MODEL>
HtlmEnsemble<MODEL>::HtlmEnsemble(const HtlmEnsembleParameters_
                                                        & params, const Geometry_ & geomTLM)
    : params_(params), ensembleSize_(params.ensembleSize),
      stateGeometry_(params.stateGeometry.value(), geomTLM.getComm()),
      incrementGeometry_(geomTLM),
      model_(stateGeometry_, params_.model.value().modelParameters.value()),
      simpleLinearModel_(geomTLM, params_.linearModel),
      controlState_(stateGeometry_, params_.control.value()),
      moderr_(stateGeometry_, params_.modelAuxControl.value()),
      modauxinc_(incrementGeometry_, params_.modelAuxIncrement.value()) {

    Log::trace() << "HtlmEnsemble<MODEL>::HtlmEnsemble() starting"
                 << std::endl;

    //  Setup perturbed member, linear Ensemble, nonlinear differences
    perturbedStates_.reserve(ensembleSize_);
    linearEnsemble_.reserve(ensembleSize_);
    nonLinearDifferences_.reserve(ensembleSize_);
    for (size_t m = 0; m < ensembleSize_; ++m) {
        perturbedStates_.emplace_back(stateGeometry_, params_.perturbedMembers.value()[m]);
        linearEnsemble_.emplace_back(incrementGeometry_, simpleLinearModel_.variables(),
                                                      controlState_.validTime());
        linearEnsemble_[m].diff(controlState_, perturbedStates_[m]);
        nonLinearDifferences_.emplace_back(incrementGeometry_, simpleLinearModel_.variables(),
                                            controlState_.validTime());
    }

    // Check you have supplied the number of members you promised
    if (params_.perturbedMembers.value().size() != ensembleSize_) {
      Log::info() << "supplied number of inital states does not equal ensemble size in yaml"
                                                                               << std::endl;
      abort();
    }

    Log::trace() << "HtlmEnsemble<MODEL>::HtlmEnsemble() done" << std::endl;
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HtlmEnsemble<MODEL>::step(const util::Duration & tstep)  {
    Log::trace() << "HtlmEnsemble<MODEL>::step() starting" << std::endl;
    State_ downsampled_Control(incrementGeometry_, controlState_);
    simpleLinearModel_.setTrajectory(controlState_, downsampled_Control, moderr_);
    PostProcessor<State_> post;
    model_.forecast(controlState_, moderr_, tstep, post);

    for (size_t m = 0; m < ensembleSize_; ++m) {
        model_.forecast(perturbedStates_[m], moderr_, tstep, post);
        simpleLinearModel_.forecastTL(linearEnsemble_[m], modauxinc_, tstep);
    }
    for (size_t m = 0; m < ensembleSize_; ++m) {
        nonLinearDifferences_[m].updateTime(tstep);
        nonLinearDifferences_[m].diff(controlState_, (perturbedStates_[m]));
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
