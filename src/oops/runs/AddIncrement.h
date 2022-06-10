/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_ADDINCREMENT_H_
#define OOPS_RUNS_ADDINCREMENT_H_

#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL>
class IncrementParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(IncrementParameters, Parameters)
  typedef Increment<MODEL> Increment_;

 public:
  typedef typename Increment_::ReadParameters_ IncrementReadParameters_;

  /// Parameters to pass to Increment::read().
  IncrementReadParameters_ read{this};

  Parameter<Variables> addedVariables{
    "added variables", "List of variables to add", {}, this};
  OptionalParameter<double> scalingFactor{
    "scaling factor", "Scaling factor for the increment", this};
};

/// YAML options taken by the AddIncrement application.
template <typename MODEL>
class AddIncrementParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(AddIncrementParameters, ApplicationParameters)
  typedef Geometry<MODEL> Geometry_;
  typedef Increment<MODEL> Increment_;
  typedef State<MODEL> State_;

 public:
  typedef typename Geometry_::Parameters_ GeometryParameters_;
  typedef IncrementParameters<MODEL> IncrementParameters_;
  typedef typename State_::Parameters_ StateParameters_;
  typedef typename State_::WriteParameters_ StateWriteParameters_;

  RequiredParameter<GeometryParameters_> stateGeometry{
    "state geometry", "State resolution", this};
  RequiredParameter<GeometryParameters_> incrementGeometry{
    "increment geometry", "Increment resolution", this};

  RequiredParameter<StateParameters_> state{
    "state", "State to be incremented", this};
  RequiredParameter<IncrementParameters_> increment{
    "increment", "Increment to add to state", this};

  RequiredParameter<StateWriteParameters_> output{
    "output", "Where to write the output", this};
};

/// Application that adds an increment to a state and writes the sum to a file.
///
/// The increment may optionally be multiplied by a scaling factor and have a different resolution
/// than the state.
template <typename MODEL> class AddIncrement : public Application {
  typedef AddIncrementParameters<MODEL> AddIncrementParameters_;
  typedef Geometry<MODEL>  Geometry_;
  typedef State<MODEL>     State_;
  typedef Increment<MODEL> Increment_;
  typedef IncrementParameters<MODEL> IncrementParameters_;

 public:
// -----------------------------------------------------------------------------
  explicit AddIncrement(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
// -----------------------------------------------------------------------------
  virtual ~AddIncrement() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
//  Load input configuration options
    AddIncrementParameters_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

//  Setup resolution
    const Geometry_ stateResol(params.stateGeometry, this->getComm());

    const Geometry_ incResol(params.incrementGeometry, this->getComm());

//  Read state
    State_ xx(stateResol, params.state);
    Log::test() << "State: " << xx << std::endl;

//  Read increment
    const IncrementParameters_ &incParams = params.increment;
    Increment_ dx(incResol, incParams.addedVariables, xx.validTime());
    dx.read(incParams.read);
    Log::test() << "Increment: " << dx << std::endl;

//  Scale increment
    if (incParams.scalingFactor.value() != boost::none) {
      dx *= *incParams.scalingFactor.value();
      Log::test() << "Scaled the increment: " << dx << std::endl;
    }

//  Assertions on state versus increment
    ASSERT(xx.validTime() == dx.validTime());

//  Add increment to state
    xx += dx;

//  Write state
    xx.write(params.output);

    Log::test() << "State plus increment: " << xx << std::endl;

    return 0;
  }
// -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    AddIncrementParameters_ params;
    params.outputSchema(outputPath);
  }
// -----------------------------------------------------------------------------
  void validateConfig(const eckit::Configuration & fullConfig) const override {
    AddIncrementParameters_ params;
    params.validate(fullConfig);
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::AddIncrement<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_ADDINCREMENT_H_
