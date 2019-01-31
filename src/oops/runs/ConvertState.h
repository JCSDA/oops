/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_CONVERTSTATE_H_
#define OOPS_RUNS_CONVERTSTATE_H_

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/State.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class ConvertState : public Application {
  typedef Geometry<MODEL> Geometry_;
  typedef State<MODEL>    State_;

 public:
// -----------------------------------------------------------------------------
  ConvertState() {}
// -----------------------------------------------------------------------------
  virtual ~ConvertState() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup resolution for intput and output
    const eckit::LocalConfiguration inputResolConfig(fullConfig, "inputresolution");
    const Geometry_ resol1(inputResolConfig);

    const eckit::LocalConfiguration outputResolConfig(fullConfig, "outputresolution");
    const Geometry_ resol2(outputResolConfig);

//  Read state
    const eckit::LocalConfiguration inputConfig(fullConfig, "input");
    std::vector<std::string> vv;
    inputConfig.get("variables", vv);
    oops::Variables vars(vv);
    State_ xx1(resol1, vars, inputConfig);
    Log::test() << "Input state: " << xx1 << std::endl;

//  Change resolution
    State_ xx2(resol2, xx1);

//  Write state
    const eckit::LocalConfiguration outputConfig(fullConfig, "output");
    xx2.write(outputConfig);

    Log::test() << "Output state: " << xx2 << std::endl;

    return 0;
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::ConvertState<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_CONVERTSTATE_H_
