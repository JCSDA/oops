/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_DIFFSTATES_H_
#define OOPS_RUNS_DIFFSTATES_H_

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the DiffStates application.
template <typename MODEL>
class DiffStatesParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(DiffStatesParameters, ApplicationParameters)

 public:
  typedef typename Geometry<MODEL>::Parameters_ GeometryParameters_;
  typedef typename State<MODEL>::Parameters_ StateParameters_;
  typedef typename Increment<MODEL>::WriteParameters_ IncrementWriteParameters_;

  /// State geometry parameters.
  RequiredParameter<GeometryParameters_> stateGeometryConf{"state geometry", this};

  /// Increment geometry parameters.
  RequiredParameter<GeometryParameters_> incGeometryConf{"increment geometry", this};

  /// First state parameters.
  RequiredParameter<StateParameters_> stateConf1{"state1", this};

  /// Second state parameters (to take away from the first).
  RequiredParameter<StateParameters_> stateConf2{"state2", this};

  /// Output increment parameters.
  RequiredParameter<IncrementWriteParameters_> outputConfig{"output", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class DiffStates : public Application {
  typedef Geometry<MODEL>  Geometry_;
  typedef State<MODEL>     State_;
  typedef Increment<MODEL> Increment_;

  typedef DiffStatesParameters<MODEL> DiffStatesParameters_;

 public:
// -----------------------------------------------------------------------------
  explicit DiffStates(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
// -----------------------------------------------------------------------------
  virtual ~DiffStates() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
//  Deserialize parameters
    DiffStatesParameters_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

//  Setup resolutions
    const Geometry_ stateGeometry(params.stateGeometryConf, this->getComm());
    const Geometry_ incGeometry(params.incGeometryConf, this->getComm());

//  Read first state
    State_ xx1(stateGeometry, params.stateConf1);
    Log::test() << "Input state 1: " << xx1 << std::endl;

//  Read second state (to take away from the first)
    State_ xx2(stateGeometry, params.stateConf2);
    Log::test() << "Input state 2: " << xx2 << std::endl;

//  Assertions on two states
    ASSERT(xx1.validTime() == xx2.validTime());

//  Create increment
    Increment_ dx(incGeometry, xx1.variables(), xx1.validTime());
    dx.diff(xx1, xx2);

//  Write increment
    dx.write(params.outputConfig);

    Log::test() << "Output increment: " << dx << std::endl;

    return 0;
  }
// -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    DiffStatesParameters_ params;
    params.outputSchema(outputPath);
  }
// -----------------------------------------------------------------------------
  void validateConfig(const eckit::Configuration & fullConfig) const override {
    DiffStatesParameters_ params;
    params.validate(fullConfig);
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::DiffStates<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_DIFFSTATES_H_
