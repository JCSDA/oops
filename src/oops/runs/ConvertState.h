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
#include "oops/base/State.h"
#include "oops/generic/instantiateVariableChangeFactory.h"
#include "oops/interface/VariableChange.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
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

  /// Variable changes
  OptionalParameter<std::vector<eckit::LocalConfiguration>> varChanges{"variable changes", this};

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
  explicit ConvertState(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateVariableChangeFactory<MODEL>();
  }
// -------------------------------------------------------------------------------------------------
  virtual ~ConvertState() {}
// -------------------------------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
//  Deserialize parameters
    ConvertStateParameters_ params;
    params.validateAndDeserialize(fullConfig);

//  Setup resolution for input and output
    const Geometry_ resol1(params.inputGeometry, this->getComm());
    const Geometry_ resol2(params.outputGeometry, this->getComm());

//  Variable transform(s)
    std::vector<VariableChange_> chvars;
    std::vector<bool> inverse;

    if (params.varChanges.value() != boost::none) {
      for (size_t cv = 0; cv < params.varChanges.value()->size(); ++cv) {
        chvars.emplace_back(resol2, (*params.varChanges.value())[cv]);
        inverse.push_back((*params.varChanges.value())[cv].getBool("do inverse", false));
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
      std::unique_ptr<State_> xx(new State_(resol2, xxi));  // Pointer that can be reset after chvar

//    Variable transform(s)
      for (size_t cv = 0; cv < chvars.size(); ++cv) {
        if (!inverse[cv]) {
          State_ xchvarout = chvars[cv].changeVar(*xx);
          xx.reset(new State_(xchvarout));
        } else {
          State_ xchvarout = chvars[cv].changeVarInverse(*xx);
          xx.reset(new State_(xchvarout));
        }
        Log::test() << "Variable transform: " << chvars[cv] << std::endl;
        Log::test() << "State after variable transform: " << *xx << std::endl;
      }

//    Write state
      xx->write(stateParams.output.value());

      Log::test() << "Output state: " << *xx << std::endl;
    }
    return 0;
  }
// -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    ConvertStateParameters_ params;
    params.outputSchema(outputPath);
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
