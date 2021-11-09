/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_GENENSPERTB_H_
#define OOPS_RUNS_GENENSPERTB_H_

#include <memory>
#include <sstream>
#include <string>


#include "eckit/config/Configuration.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/Model.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/base/StateWriter.h"
#include "oops/base/Variables.h"
#include "oops/interface/ModelAuxControl.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

/// Options taken by the GenEnsPertB application.
template <typename MODEL> class GenEnsPertBParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(GenEnsPertBParameters, ApplicationParameters)

 public:
  typedef ModelSpaceCovarianceParametersWrapper<MODEL> CovarianceParameters_;
  typedef typename Geometry<MODEL>::Parameters_        GeometryParameters_;
  typedef ModelParametersWrapper<MODEL>                ModelParameters_;
  typedef State<MODEL>                                 State_;
  typedef StateParametersND<MODEL>                     StateParametersND_;
  typedef StateWriterParameters<State_>                StateWriterParameters_;

  /// Geometry parameters.
  RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Model parameters.
  RequiredParameter<ModelParameters_> model{"model", this};

  /// Initial state parameters.
  RequiredParameter<StateParametersND_> initialCondition{"initial condition", this};

  /// Augmented model state.
  Parameter<eckit::LocalConfiguration> modelAuxControl{
    "model aux control", eckit::LocalConfiguration(), this};

  /// Forecast length.
  RequiredParameter<util::Duration> forecastLength{"forecast length", this};

  /// List of variables to perturb.
  RequiredParameter<Variables> perturbedVariables{"perturbed variables", this};

  /// Background error covariance model.
  RequiredParameter<CovarianceParameters_> backgroundError{"background error", this};

  /// Size of the perturbed ensemble to generate.
  RequiredParameter<int> members{"members", this};

  /// Where to write the output.
  RequiredParameter<StateWriterParameters_> output{"output", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class GenEnsPertB : public Application {
  typedef ModelSpaceCovarianceBase<MODEL>           CovarianceBase_;
  typedef CovarianceFactory<MODEL>                  CovarianceFactory_;
  typedef ModelSpaceCovarianceParametersBase<MODEL> CovarianceParametersBase_;
  typedef Geometry<MODEL>                           Geometry_;
  typedef Model<MODEL>                              Model_;
  typedef ModelAuxControl<MODEL>                    ModelAux_;
  typedef Increment<MODEL>                          Increment_;
  typedef State<MODEL>                              State_;
  typedef StateWriterParameters<State_>             StateWriterParameters_;

  typedef GenEnsPertBParameters<MODEL>              GenEnsPertBParameters_;

 public:
// -----------------------------------------------------------------------------
  explicit GenEnsPertB(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateCovarFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~GenEnsPertB() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
//  Deserialize parameters
    GenEnsPertBParameters_ params;
    params.validateAndDeserialize(fullConfig);

//  Setup resolution
    const Geometry_ resol(params.geometry, this->getComm(), oops::mpi::myself());

//  Setup Model
    const Model_ model(resol, params.model.value().modelParameters);

//  Setup initial state
    const State_ xx(resol, params.initialCondition);
    Log::test() << "Initial state: " << xx << std::endl;

//  Setup augmented state
    const ModelAux_ moderr(resol, params.modelAuxControl);

//  Setup times
    const util::Duration fclength = params.forecastLength;
    const util::DateTime bgndate(xx.validTime());
    const util::DateTime enddate(bgndate + fclength);
    Log::info() << "Running forecast from " << bgndate << " to " << enddate << std::endl;

//  Setup variables
    const Variables &vars = params.perturbedVariables;

//  Setup B matrix
    const CovarianceParametersBase_ &covarParams =
        params.backgroundError.value().covarianceParameters;
    std::unique_ptr<CovarianceBase_> Bmat(CovarianceFactory_::create(
                                            covarParams, resol, vars, xx, xx));

//  Generate perturbed states
    Increment_ dx(resol, vars, bgndate);
    for (int jm = 0; jm < params.members; ++jm) {
//    Generate pertubation
      Bmat->randomize(dx);

//    Add mean state
      State_ xp(xx);
      xp += dx;

//    Setup forecast outputs
      PostProcessor<State_> post;

      StateWriterParameters_ outParams = params.output;
      outParams.write.setMember(jm + 1);

      post.enrollProcessor(new StateWriter<State_>(outParams));

//    Run forecast
      model.forecast(xp, moderr, fclength, post);
      Log::test() << "Member " << jm << " final state: " << xp << std::endl;
    }

    return 0;
  }
// -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    GenEnsPertBParameters_ params;
    params.outputSchema(outputPath);
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::GenEnsPertB<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_GENENSPERTB_H_
