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

#include "eckit/config/LocalConfiguration.h"
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

  // Constructor declared for convenience; it can be removed once top-level parameters
  // for the AddIncrement application are defined.
  IncrementParameters(const eckit::Configuration &conf, const std::string &path) {
    validateAndDeserialize(eckit::LocalConfiguration(conf, path));
  }

  /// List of variables to add.
  Parameter<Variables> addedVariables{"added variables", {}, this};
  /// Parameters to pass to Increment::read().
  IncrementReadParameters_ read{this};
  /// Scaling factor for the increment.
  OptionalParameter<double> scalingFactor{"scaling factor", this};
};

template <typename MODEL> class AddIncrement : public Application {
  typedef Geometry<MODEL>  Geometry_;
  typedef State<MODEL>     State_;
  typedef Increment<MODEL> Increment_;

 public:
// -----------------------------------------------------------------------------
  explicit AddIncrement(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
// -----------------------------------------------------------------------------
  virtual ~AddIncrement() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup resolution
    const eckit::LocalConfiguration stateResolConf(fullConfig, "state geometry");
    const Geometry_ stateResol(stateResolConf, this->getComm());

    const eckit::LocalConfiguration incResolConf(fullConfig, "increment geometry");
    const Geometry_ incResol(incResolConf, this->getComm());

//  Read state
    const eckit::LocalConfiguration stateConf(fullConfig, "state");
    State_ xx(stateResol, stateConf);
    Log::test() << "State: " << xx << std::endl;

//  Read increment
    IncrementParameters<MODEL> incParams(fullConfig, "increment");
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
    const eckit::LocalConfiguration outputConfig(fullConfig, "output");
    xx.write(outputConfig);

    Log::test() << "State plus increment: " << xx << std::endl;

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::AddIncrement<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_ADDINCREMENT_H_
