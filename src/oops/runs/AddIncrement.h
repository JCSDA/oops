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
#include "oops/base/State.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

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
    const eckit::LocalConfiguration incConf(fullConfig, "increment");
    std::vector<std::string> incvv;
    incConf.get("added variables", incvv);
    oops::Variables incVars(incvv);
    Increment_ dx(incResol, incVars, xx.validTime());
    dx.read(incConf);
    Log::test() << "Increment: " << dx << std::endl;

//  Scale increment
    if (incConf.has("scaling factor")) {
      dx *= incConf.getDouble("scaling factor");
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
