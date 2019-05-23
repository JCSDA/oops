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
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class DiffStates : public Application {
  typedef Geometry<MODEL>  Geometry_;
  typedef State<MODEL>     State_;
  typedef Increment<MODEL> Increment_;

 public:
// -----------------------------------------------------------------------------
  DiffStates() {}
// -----------------------------------------------------------------------------
  virtual ~DiffStates() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup resolution
    const eckit::LocalConfiguration stateResolConf(fullConfig, "stateResol");
    const Geometry_ stateResol(stateResolConf);

    const eckit::LocalConfiguration incResolConf(fullConfig, "incrementResol");
    const Geometry_ incResol(incResolConf);

//  Read first state
    const eckit::LocalConfiguration stateConf1(fullConfig, "state1");
    std::vector<std::string> vv1;
    stateConf1.get("variables", vv1);
    oops::Variables vars1(vv1);
    State_ xx1(stateResol, vars1, stateConf1);
    Log::test() << "Input state 1: " << xx1 << std::endl;

//  Read second state (to take away from the first)
    const eckit::LocalConfiguration stateConf2(fullConfig, "state2");
    std::vector<std::string> vv2;
    stateConf2.get("variables", vv2);
    oops::Variables vars2(vv2);
    State_ xx2(stateResol, vars2, stateConf2);
    Log::test() << "Input state 2: " << xx2 << std::endl;

//  Assertions on two states
    ASSERT(xx1.validTime() == xx2.validTime());

//  Create increment
    Increment_ dx(incResol, vars1, xx1.validTime());
    dx.diff(xx1, xx2);

//  Write increment
    const eckit::LocalConfiguration outputConfig(fullConfig, "output");
    dx.write(outputConfig);

    Log::test() << "Output increment: " << dx << std::endl;

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::DiffStates<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_DIFFSTATES_H_
