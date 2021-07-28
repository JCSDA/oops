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

namespace oops {

template <typename MODEL> class ConvertState : public Application {
  typedef Geometry<MODEL>              Geometry_;
  typedef State<MODEL>                 State_;
  typedef VariableChange<MODEL>    VariableChange_;

 public:
// -------------------------------------------------------------------------------------------------
  explicit ConvertState(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateVariableChangeFactory<MODEL>();
  }
// -------------------------------------------------------------------------------------------------
  virtual ~ConvertState() {}
// -------------------------------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup resolution for intput and output
    const eckit::LocalConfiguration inputResolConfig(fullConfig, "input geometry");
    const Geometry_ resol1(inputResolConfig, this->getComm());

    const eckit::LocalConfiguration outputResolConfig(fullConfig, "output geometry");
    const Geometry_ resol2(outputResolConfig, this->getComm());

//  Variable transform(s)
    std::vector<VariableChange_> chvars;
    std::vector<bool> inverse;

    std::vector<eckit::LocalConfiguration> chvarconfs;
    fullConfig.get("variable changes", chvarconfs);
    for (size_t cv = 0; cv < chvarconfs.size(); ++cv) {
      chvars.emplace_back(resol2, chvarconfs[cv]);
      inverse.push_back(chvarconfs[cv].getBool("do inverse", false));
    }

//  List of input and output states
    std::vector<eckit::LocalConfiguration> statesConf;
    fullConfig.get("states", statesConf);
    int nstates = statesConf.size();

//  Loop over states
    for (int jm = 0; jm < nstates; ++jm) {
//    Print output
      Log::info() << "Converting state " << jm+1 << " of " << nstates << std::endl;

//    Read state
      const eckit::LocalConfiguration inputConfig(statesConf[jm], "input");
      State_ xxi(resol1, inputConfig);
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
      const eckit::LocalConfiguration outputConfig(statesConf[jm], "output");
      xx->write(outputConfig);

      Log::test() << "Output state: " << *xx << std::endl;
    }
    return 0;
  }
// -------------------------------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::ConvertState<" + MODEL::name() + ">";
  }
// -------------------------------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_CONVERTSTATE_H_
