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
#include "oops/interface/VariableChange.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/Logger.h"

namespace oops {

// -------------------------------------------------------------------------------------------------

template <typename MODEL> class ConvertState : public Application {
  typedef Geometry<MODEL>               Geometry_;
  typedef State<MODEL>                  State_;
  typedef VariableChange<MODEL>         VariableChange_;

 public:
// -------------------------------------------------------------------------------------------------
  explicit ConvertState(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
// -------------------------------------------------------------------------------------------------
  virtual ~ConvertState() {}
// -------------------------------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
//  Setup resolution for input and output
    const Geometry_ resol1(eckit::LocalConfiguration(fullConfig, "input geometry"),
                           this->getComm());
    const Geometry_ resol2(eckit::LocalConfiguration(fullConfig, "output geometry"),
                           this->getComm());

    // Setup change of variable
    std::unique_ptr<VariableChange_> vc;
    oops::Variables varout;
    bool inverse = false;
    const eckit::LocalConfiguration chconf = fullConfig.getSubConfiguration("variable change");
    if (chconf.has("output variables")) {
      vc.reset(new VariableChange_(chconf, resol2));
      varout = Variables(chconf, "output variables");
      inverse = chconf.getBool("do inverse", false);
    }

//  List of input and output states
    const std::vector<eckit::LocalConfiguration> stateConfs =
        fullConfig.getSubConfigurations("states");
    const int nstates = stateConfs.size();

//  Loop over states
    for (int jm = 0; jm < nstates; ++jm) {
//    Print output
      Log::info() << "Converting state " << jm+1 << " of " << nstates << std::endl;

//    Read state
      State_ xxi(resol1, eckit::LocalConfiguration(stateConfs[jm], "input"));
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
      eckit::LocalConfiguration outconf(stateConfs[jm], "output");
      xx.write(outconf);

      Log::test() << "Output state: " << xx << std::endl;
    }
    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::ConvertState<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_CONVERTSTATE_H_
