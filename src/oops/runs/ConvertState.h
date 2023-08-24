/*
 * (C) Copyright 2018-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_CONVERTSTATE_H_
#define OOPS_RUNS_CONVERTSTATE_H_

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Geometry.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/State.h"
#include "oops/interface/VariableChange.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

template <typename MODEL> class ConvertStateStatesParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(ConvertStateStatesParameters, Parameters)
  typedef State<MODEL> State_;

 public:
  typedef typename State_::Parameters_      StateParameters_;
  typedef typename State_::WriteParameters_ WriteParameters_;

  RequiredParameter<StateParameters_> input{"input", this};
  RequiredParameter<WriteParameters_> output{"output", this};
};

/// Options controlling variable change in the ConvertState application.
template <typename MODEL> class VarChangeParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(VarChangeParameters, Parameters)
  typedef typename VariableChange<MODEL>::Parameters_ VariableChangeParameters_;

 public:
  // parameters for variable change.
  VariableChangeParameters_ varChange{this};
  Parameter<bool> doInverse{"do inverse",
                            "apply inverse variable change instead of variable change",
                            false, this};
};

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

/// Options taken by the ConvertState application.
template <typename MODEL> class ConvertStateParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(ConvertStateParameters, ApplicationParameters)
  typedef Geometry<MODEL> Geometry_;

 public:
  typedef typename Geometry_::Parameters_ GeometryParameters_;

  /// Input Geometry parameters.
  RequiredParameter<GeometryParameters_> inputGeometry{"input geometry", this};

  /// Output Geometry parameters.
  RequiredParameter<GeometryParameters_> outputGeometry{"output geometry", this};

  /// Variable change parameters (and option to do inverse).
  OptionalParameter<VarChangeParameters<MODEL>> varChange{"variable change", this};

  /// States to be converted
  RequiredParameter<std::vector<ConvertStateStatesParameters<MODEL>>> states{"states", this};
};

// -------------------------------------------------------------------------------------------------
// -------------------------------------------------------------------------------------------------

template <typename MODEL> class ConvertState : public Application {
  typedef Geometry<MODEL>               Geometry_;
  typedef State<MODEL>                  State_;
  typedef VariableChange<MODEL>         VariableChange_;
  typedef ConvertStateParameters<MODEL> ConvertStateParameters_;
  typedef ConvertStateStatesParameters<MODEL> ConvertStateStatesParameters_;

 public:
// -------------------------------------------------------------------------------------------------
  explicit ConvertState(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
// -------------------------------------------------------------------------------------------------
  virtual ~ConvertState() {}
// -------------------------------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
//  Deserialize parameters
    ConvertStateParameters_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

//  Setup resolution for input and output
    const Geometry_ resol1(params.inputGeometry, this->getComm());
    const Geometry_ resol2(params.outputGeometry, this->getComm());

    // Setup change of variable
    std::unique_ptr<VariableChange_> vc;
    oops::Variables varout;
    bool inverse = false;
    if (params.varChange.value() != boost::none) {
      eckit::LocalConfiguration chconf(params.varChange.value()->toConfiguration());
      if (chconf.has("output variables")) {
        vc.reset(new VariableChange_(chconf, resol2));
        varout = Variables(chconf, "output variables");
        inverse = chconf.getBool("do inverse", false);
      }
    }

//  List of input and output states
    const int nstates = params.states.value().size();

//  Loop over states
    for (int jm = 0; jm < nstates; ++jm) {
//    Read current state parameters
      const ConvertStateStatesParameters_ stateParams = params.states.value()[jm];

//    Print output
      Log::info() << "Converting state " << jm+1 << " of " << nstates << std::endl;

//    Read state
      State_ xxi(resol1, stateParams.input.value());
      Log::test() << "Input state: " << xxi << std::endl;

//    Copy and change resolution
      State_ xx(resol2, xxi);

//    Variable transform(s)
      if (vc) {
          // Create variable change
        oops::Variables varin = xx.variables();
        if (inverse) {
          vc->changeVarInverse(xx, varout);
        } else {
          vc->changeVar(xx, varout);
        }
        Log::test() << "Variable transform: " << *vc << std::endl;
        Log::test() << "Variable change from: " << varin << std::endl;
        Log::test() << "Variable change to: " << varout << std::endl;
        Log::test() << "State after variable transform: " << xx << std::endl;
      }

//    Write state
      xx.write(stateParams.output.value());

      Log::test() << "Output state: " << xx << std::endl;
    }
    return 0;
  }
// -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    ConvertStateParameters_ params;
    params.outputSchema(outputPath);
  }
// -----------------------------------------------------------------------------
  void validateConfig(const eckit::Configuration & fullConfig) const override {
    ConvertStateParameters_ params;
    params.validate(fullConfig);
  }
// -------------------------------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::ConvertState<" + MODEL::name() + ">";
  }
// -------------------------------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_CONVERTSTATE_H_
