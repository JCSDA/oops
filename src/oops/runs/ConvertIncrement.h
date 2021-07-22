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
#include "oops/base/LinearVariableChangeBase.h"
#include "oops/generic/instantiateVariableChangeFactory.h"
#include "oops/interface/State.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class ConvertIncrement : public Application {
  typedef Geometry<MODEL>                    Geometry_;
  typedef Increment<MODEL>                   Increment_;
  typedef State<MODEL>                       State_;
  typedef LinearVariableChangeBase<MODEL>    LinearVariableChange_;
  typedef LinearVariableChangeFactory<MODEL> LinearVariableChangeFactory_;

 public:
// -------------------------------------------------------------------------------------------------
  explicit ConvertIncrement(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm)
  {
    instantiateVariableChangeFactory<MODEL>();
  }
// -------------------------------------------------------------------------------------------------
  virtual ~ConvertIncrement() {}
// -------------------------------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const {
//  Setup resolution for intput and output
    const eckit::LocalConfiguration inputResolConfig(fullConfig, "input geometry");
    const Geometry_ resol1(inputResolConfig, this->getComm());

    const eckit::LocalConfiguration outputResolConfig(fullConfig, "output geometry");
    const Geometry_ resol2(outputResolConfig, this->getComm());

//  Variable transform(s)
    std::vector<bool> inverse;
    std::vector<bool> adjoint;

    std::vector<eckit::LocalConfiguration> chvarconfs;
    fullConfig.get("linear variable changes", chvarconfs);
    for (size_t cv = 0; cv < chvarconfs.size(); ++cv) {
      inverse.push_back(chvarconfs[cv].getBool("do inverse", false));
      adjoint.push_back(chvarconfs[cv].getBool("do adjoint", false));
    }

//  List of input and output increments
    std::vector<eckit::LocalConfiguration> incrementsConf;
    fullConfig.get("increments", incrementsConf);
    int nincrements = incrementsConf.size();

//  Loop over increments
    for (int jm = 0; jm < nincrements; ++jm) {
//    Print output
      Log::info() << "Converting increment " << jm+1 << " of " << nincrements << std::endl;

//    Datetime for incrmement
      const util::DateTime incdatetime(incrementsConf[jm].getString("date"));

//    Variables for input increment
      const Variables incvars(incrementsConf[jm], "input variables");

//    Read input
      const eckit::LocalConfiguration inputConfig(incrementsConf[jm], "input");
      Increment_ dxi(resol1, incvars, incdatetime);
      dxi.read(inputConfig);
      Log::test() << "Input increment: " << dxi << std::endl;

//    Copy and change resolution
      std::unique_ptr<Increment_> dx(new Increment_(resol2, dxi));  // Pointer that can be reset

//    Trajectory state for linear variable transform
      std::unique_ptr<State_> xtraj;  // Pointer that can be reset

//    Variable transform(s)
      for (size_t cv = 0; cv < chvarconfs.size(); ++cv) {
        // Read trajectory
        if (cv == 0) {
          const eckit::LocalConfiguration trajConfig(incrementsConf[jm], "trajectory");
          xtraj.reset(new State_(resol1, trajConfig));
          ASSERT(xtraj->validTime() == dx->validTime());  // Check time is consistent
          Log::test() << "Trajectory state: " << *xtraj << std::endl;
        }

        // Create variable change
        std::unique_ptr<LinearVariableChange_> lvc;
        lvc.reset(LinearVariableChangeFactory_::create(*xtraj, *xtraj, resol2, chvarconfs[cv]));

        // Print info
        Log::info() << "Variable transform " << cv+1 << " of " << chvarconfs.size() << ": "
                    << *lvc << std::endl;

        Increment_ xchvarout = lvc->multiply(*dx);
        dx.reset(new Increment_(xchvarout));

        Log::test() << "Increment after variable transform: " << *dx << std::endl;
      }

//    Write state
      const eckit::LocalConfiguration outputConfig(incrementsConf[jm], "output");
      dx->write(outputConfig);

      Log::test() << "Output increment: " << *dx << std::endl;
    }
    return 0;
  }
// -------------------------------------------------------------------------------------------------
 private:
  std::string appname() const {
    return "oops::ConvertIncrement<" + MODEL::name() + ">";
  }
// -------------------------------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_CONVERTINCREMENT_H_
