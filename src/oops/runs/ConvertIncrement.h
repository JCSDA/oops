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
#include "oops/base/State.h"
#include "oops/generic/instantiateVariableChangeFactory.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

/// Options describing each increment read by the ConvertIncrement application.
template <typename MODEL> class IncrementParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(IncrementParameters, Parameters);

 public:
  typedef typename Increment<MODEL>::ReadParameters_  ReadParameters_;
  typedef typename Increment<MODEL>::WriteParameters_ WriteParameters_;
  typedef typename State<MODEL>::Parameters_          StateParameters_;

  RequiredParameter<util::DateTime> date{"date", this};
  RequiredParameter<oops::Variables> inputVariables{"input variables", this};
  RequiredParameter<ReadParameters_> input{"input", this};
  RequiredParameter<WriteParameters_> output{"output", this};
  RequiredParameter<StateParameters_> trajectory{"trajectory", this};
};

/// Top-level options taken by the ConvertIncrement application.
template <typename MODEL> class ConvertIncrementParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(ConvertIncrementParameters, ApplicationParameters);

 public:
  typedef typename Geometry<MODEL>::Parameters_ GeometryParameters_;
  typedef IncrementParameters<MODEL>            IncrementParameters_;

  /// Input geometry parameters.
  RequiredParameter<GeometryParameters_> inputGeometry{"input geometry", this};

  /// Output geometry parameters.
  RequiredParameter<GeometryParameters_> outputGeometry{"output geometry", this};

  /// List of linear variable changes.
  RequiredParameter<std::vector<eckit::LocalConfiguration>> linearVarChanges{
      "linear variable changes", this};

  /// List of increments.
  RequiredParameter<std::vector<IncrementParameters_>> increments{"increments", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class ConvertIncrement : public Application {
  typedef Geometry<MODEL>                    Geometry_;
  typedef Increment<MODEL>                   Increment_;
  typedef State<MODEL>                       State_;
  typedef LinearVariableChangeBase<MODEL>    LinearVariableChange_;
  typedef LinearVariableChangeFactory<MODEL> LinearVariableChangeFactory_;

  typedef typename Increment<MODEL>::ReadParameters_  ReadParameters_;
  typedef typename Increment<MODEL>::WriteParameters_ WriteParameters_;
  typedef IncrementParameters<MODEL>                  IncrementParameters_;

  typedef ConvertIncrementParameters<MODEL>  Parameters_;

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
//  Deserialize parameters
    Parameters_ params;
    params.validateAndDeserialize(fullConfig);

//  Setup resolution for intput and output
    const Geometry_ resol1(params.inputGeometry, this->getComm());
    const Geometry_ resol2(params.outputGeometry, this->getComm());

//  Variable transform(s)
    std::vector<bool> inverse;
    std::vector<bool> adjoint;

    const std::vector<eckit::LocalConfiguration>& chvarconfs = params.linearVarChanges;
    for (size_t cv = 0; cv < chvarconfs.size(); ++cv) {
      inverse.push_back(chvarconfs[cv].getBool("do inverse", false));
      adjoint.push_back(chvarconfs[cv].getBool("do adjoint", false));
    }

//  List of input and output increments
    const std::vector<IncrementParameters_>& incrementParams = params.increments;
    const int nincrements = incrementParams.size();

//  Loop over increments
    for (int jm = 0; jm < nincrements; ++jm) {
//    Print output
      Log::info() << "Converting increment " << jm+1 << " of " << nincrements << std::endl;

//    Datetime for increment
      const util::DateTime incdatetime = incrementParams[jm].date;

//    Variables for input increment
      const Variables incvars = incrementParams[jm].inputVariables;

//    Read input
      const ReadParameters_ inputParams = incrementParams[jm].input;
      Increment_ dxi(resol1, incvars, incdatetime);
      dxi.read(inputParams);
      Log::test() << "Input increment: " << dxi << std::endl;

//    Copy and change resolution
      std::unique_ptr<Increment_> dx(new Increment_(resol2, dxi));  // Pointer that can be reset

//    Trajectory state for linear variable transform
      std::unique_ptr<State_> xtraj;  // Pointer that can be reset

//    Variable transform(s)
      for (size_t cv = 0; cv < chvarconfs.size(); ++cv) {
        // Read trajectory
        if (cv == 0) {
          xtraj.reset(new State_(resol1, incrementParams[jm].trajectory));
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
      const WriteParameters_ outputParams = incrementParams[jm].output;
      dx->write(outputParams);

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
