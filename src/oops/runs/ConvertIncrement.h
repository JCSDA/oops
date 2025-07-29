/*
 * (C) Copyright 2018-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_CONVERTINCREMENT_H_
#define OOPS_RUNS_CONVERTINCREMENT_H_

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/interface/LinearVariableChange.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL> class ConvertIncrement : public Application {
  typedef Geometry<MODEL>                    Geometry_;
  typedef Increment<MODEL>                   Increment_;
  typedef State<MODEL>                       State_;
  typedef LinearVariableChange<MODEL>        LinearVariableChange_;

 public:
// -------------------------------------------------------------------------------------------------
  explicit ConvertIncrement(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {}
// -------------------------------------------------------------------------------------------------
  virtual ~ConvertIncrement() {}
// -------------------------------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
//  Setup resolution for intput and output
    const Geometry_ resol1(eckit::LocalConfiguration(fullConfig, "input geometry"),
                           this->getComm());
    const Geometry_ resol2(eckit::LocalConfiguration(fullConfig, "output geometry"),
                           this->getComm());

//  Check if there is a change of variable defined in Parameters
    const eckit::LocalConfiguration linCVconf =
        fullConfig.getSubConfiguration("linear variable change");
    bool lvcDefined = linCVconf.has("output variables");

//  List of input and output increments
    const std::vector<eckit::LocalConfiguration> incConfs =
        fullConfig.getSubConfigurations("increments");
    const int nincrements = incConfs.size();

//  Loop over increments
    for (int jm = 0; jm < nincrements; ++jm) {
//    Print output
      Log::info() << "Converting increment " << jm+1 << " of " << nincrements << std::endl;

//    Datetime for increment
      const util::DateTime incdatetime(incConfs[jm].getString("date"));

//    Variables for input increment
      const Variables incvars(incConfs[jm], "input variables");

//    Read input
      const eckit::LocalConfiguration inputParams(incConfs[jm], "input");
      Increment_ dxi(resol1, incvars, incdatetime);
      dxi.read(inputParams);
      Log::test() << "Input increment: " << dxi << std::endl;

//    Copy and change resolution
      Increment_ dx(resol2, dxi);

//    Variable transform
      if (lvcDefined) {
        const eckit::LocalConfiguration trajConf(incConfs[jm], "trajectory");
        State_ xTrajBg(resol1, trajConf);
        ASSERT(xTrajBg.validTime() == dx.validTime());  // Check time is consistent
        Log::test() << "Trajectory state: " << xTrajBg << std::endl;

        // Create variable change
        LinearVariableChange_ lvc(resol2, linCVconf);
        Variables varout(linCVconf, "output variables");
        lvc.changeVarTraj(xTrajBg, varout);
        if (linCVconf.getBool("do inverse", false)) {
          lvc.changeVarInverseTL(dx, varout);
        } else {
          lvc.changeVarTL(dx, varout);
        }
      }

//    Write state
      const eckit::LocalConfiguration outputParams(incConfs[jm], "output");
      dx.write(outputParams);

      Log::test() << "Output increment: " << dx << std::endl;
    }
    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::ConvertIncrement<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_CONVERTINCREMENT_H_
