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
#include "oops/interface/ModelAuxControl.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Expect.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"


namespace oops {

/// paramaters for for forwarding the ensemble.
template <typename MODEL>
class HtlmEnsembleParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(HtlmEnsembleParameters, Parameters);
  typedef typename State<MODEL>::Parameters_    StateParameters_;
  typedef ModelParametersWrapper<MODEL>         ModelParameters_;
  typedef typename Geometry<MODEL>::Parameters_ GeometryParameters_;

 public:
  /// Nonlinear forecast model.
  RequiredParameter<ModelParameters_> model{"model", this};
  /// Background and analysis geometry.
  RequiredParameter<GeometryParameters_> stateGeometry{"state geometry", this};
  /// Configurations of all perturbed members.
  RequiredParameter<std::vector<StateParameters_>>
        perturbedMembers{"perturbed members", this};
  /// Configurations of the control member.
  RequiredParameter<StateParameters_> control{"control member", this};
  /// Augmented model state
  Parameter<eckit::LocalConfiguration> modelAuxControl{"model aux control",
    eckit::LocalConfiguration(), this};
  /// Number of perturbed ensemble members
  RequiredParameter<size_t> ensembleSize{"ensemble size", this};
};

// class declerations for htlmensemble
template <typename MODEL>
class HtlmEnsemble{
  typedef Geometry<MODEL>                               Geometry_;
  typedef Model<MODEL>                                  Model_;
  typedef ModelAuxControl<MODEL>                        ModelAux_;
  typedef State<MODEL>                                  State_;
  typedef Increment<MODEL>                              Increment_;
  typedef HtlmEnsembleParameters<MODEL>                 HtlmEnsembleParameters_;

 public:
  static const std::string classname() {return "oops::HtlmEnsemble";}

/// constructor
    HtlmEnsemble(const HtlmEnsembleParameters_ &, const Geometry_ &);

/// step forward and get coeffs
    void step(const util::Duration &) const;

 private:
  const HtlmEnsembleParameters_ params_;
  const size_t ensembleSize_;
  const Geometry_ stateGeometry_;

///  Model
  const Model_ model_;

///  control member IC
  std::unique_ptr<State_> controlState_;
///  perturbed member ICs
  std::vector<std::unique_ptr<State_>> perturbedStates_;
///  Augmented state
  std::unique_ptr<ModelAux_> moderr_;
};

//----------------------------------------------------------------------------------

template<typename MODEL>
HtlmEnsemble<MODEL>::HtlmEnsemble(const HtlmEnsembleParameters_
                                                        & params, const Geometry_ & geom)
    : params_(params), ensembleSize_(params.ensembleSize),
      stateGeometry_(params.stateGeometry.value(), geom.getComm()),
      model_(stateGeometry_, params_.model.value().modelParameters.value()) {

    Log::trace() << "HtlmEnsemble<MODEL>::HtlmEnsemble() starting"
                 << std::endl;

    ///  Setup control member IC
    controlState_ = std::make_unique<State_>(stateGeometry_, params_.control.value());
    ///  Setup perturbed member ICs
    perturbedStates_.reserve(ensembleSize_);
    for (size_t m = 0; m < ensembleSize_; ++m) {
        perturbedStates_.emplace_back(std::make_unique<State_>(stateGeometry_,
                                                              params_.perturbedMembers.
                                                              value()[m]));
    }

     /// Check you have supplied the number of members you promised
    if (params_.perturbedMembers.value().size() != ensembleSize_) {
      Log::info() << "supplied number of inital states does not equal ensemble size in yaml"
                                                                               << std::endl;
      abort();
    }

    ///  Setup augmented state
    moderr_ = std::make_unique<ModelAux_>(stateGeometry_, params_.modelAuxControl.value());
    Log::trace() << "HtlmEnsemble<MODEL>::HtlmEnsemble() done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void HtlmEnsemble<MODEL>::step(const util::Duration & tstep) const {
    Log::trace() << "HtlmEnsemble<MODEL>::step() starting" << std::endl;
    PostProcessor<State_> post;
    model_.forecast((*controlState_), *moderr_, tstep, post);
    for (size_t m = 0; m < ensembleSize_; ++m) {
        model_.forecast(*(perturbedStates_[m]), *moderr_, tstep, post);
        }

    Log::trace() << "HtlmEnsemble<MODEL>::step() done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_HTLMENSEMBLE_H_
