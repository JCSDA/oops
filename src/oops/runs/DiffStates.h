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

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL> class DiffStates : public Application {
  typedef Geometry<MODEL>  Geometry_;
  typedef State<MODEL>     State_;
  typedef Increment<MODEL> Increment_;

 public:
// -----------------------------------------------------------------------------
  explicit DiffStates(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
// -----------------------------------------------------------------------------
  virtual ~DiffStates() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
//  Setup resolutions
    const Geometry_ stateGeometry(eckit::LocalConfiguration(fullConfig, "state geometry"),
                                  this->getComm());
    const Geometry_ incGeometry(eckit::LocalConfiguration(fullConfig, "increment geometry"),
                                this->getComm());

//  Read first state
    State_ xx1(stateGeometry, eckit::LocalConfiguration(fullConfig, "state1"));
    Log::test() << "Input state 1: " << xx1 << std::endl;

//  Read second state (to take away from the first)
    State_ xx2(stateGeometry, eckit::LocalConfiguration(fullConfig, "state2"));
    Log::test() << "Input state 2: " << xx2 << std::endl;

//  Assertions on two states
    ASSERT(xx1.validTime() == xx2.validTime());

//  Create increment
    Increment_ dx(incGeometry, xx1.variables(), xx1.validTime());
    dx.diff(xx1, xx2);

//  Write increment
    dx.write(eckit::LocalConfiguration(fullConfig, "output"));

    Log::test() << "Output increment: " << dx << std::endl;

    return 0;
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
