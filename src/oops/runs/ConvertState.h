/*
 * (C) Copyright 2018 UCAR
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
#include "oops/base/VariableChangeBase.h"
#include "oops/generic/instantiateVariableChangeFactory.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/State.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class ConvertState : public Application {
  typedef Geometry<MODEL>              Geometry_;
  typedef State<MODEL>                 State_;
  typedef VariableChangeBase<MODEL>    VariableChange_;
  typedef VariableChangeFactory<MODEL> VariableChangeFactory_;

 public:
// -------------------------------------------------------------------------------------------------
  explicit ConvertState(const eckit::mpi::Comm & comm = oops::mpi::comm()) : Application(comm) {
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

//  Output variables
    oops::Variables varsout(fullConfig, "output variables");

//  Variable transform, identity if not specified in config
    std::unique_ptr<VariableChange_> changevar(VariableChangeFactory_::create(fullConfig, resol2));
    bool inverse = fullConfig.getBool("do inverse", false);

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
      State_ xx(resol2, xxi);

//    New state with variables after variable change
      State_ xxo(resol2, varsout, xxi.validTime());

//    Variable transform
      if (!inverse) {
        changevar->changeVar(xx, xxo);
      } else {
        changevar->changeVarInverse(xx, xxo);
      }

//    Write state
      const eckit::LocalConfiguration outputConfig(statesConf[jm], "output");
      xxo.write(outputConfig);

      Log::test() << "Output state: " << xxo << std::endl;
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
