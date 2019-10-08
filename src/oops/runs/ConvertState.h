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
    const eckit::LocalConfiguration inputResolConfig(fullConfig, "inputresolution");
    const Geometry_ resol1(inputResolConfig, this->getComm());

    const eckit::LocalConfiguration outputResolConfig(fullConfig, "outputresolution");
    const Geometry_ resol2(outputResolConfig, this->getComm());

//  Read state
    const eckit::LocalConfiguration inputConfig(fullConfig, "input");
    const std::vector<std::string> vi = inputConfig.getStringVector("variables");
    oops::Variables vars(vi);
    State_ xx1(resol1, vars, inputConfig);
    Log::test() << "Input state: " << xx1 << std::endl;

//  Change resolution
    State_ xx2(resol2, xx1);

//  New variables if any
    const std::vector<std::string> vo = fullConfig.getStringVector("outputVariables.variables", vi);
    oops::Variables varsnew(vo);

//  New state with variables after variable change
    State_ xx3(resol2, varsnew, xx2.validTime());

//  Perform transform
    std::unique_ptr<VariableChange_> changevar(VariableChangeFactory_::create(fullConfig, resol2));
    bool inverse = fullConfig.getBool("doinverse", false);
    if (!inverse) {
      changevar->changeVar(xx2, xx3);
    } else {
      changevar->changeVarInverse(xx2, xx3);
    }

//  Write state
    const eckit::LocalConfiguration outputConfig(fullConfig, "output");
    xx3.write(outputConfig);

    Log::test() << "Output state: " << xx3 << std::endl;

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
