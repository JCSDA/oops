/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_CONVERTTOSTRUCTUREDGRID_H_
#define OOPS_RUNS_CONVERTTOSTRUCTUREDGRID_H_

#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/State.h"
#include "oops/base/StateEnsemble.h"
#include "oops/base/StructuredGridWriter.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL> class StateToStructuredGridParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(StateToStructuredGridParameters, ApplicationParameters)
  typedef State<MODEL>                   State_;
  typedef Geometry<MODEL>                Geometry_;

 public:
  typedef typename Geometry_::Parameters_   GeometryParameters_;

  RequiredParameter<GeometryParameters_>       stateGeometry{"state geometry", this};
  RequiredParameter<eckit::LocalConfiguration> state{"state", this};
  RequiredParameter<eckit::LocalConfiguration> structuredGridInterp
                   {"structured grid interpolation", this};
};

template <typename MODEL> class IncToStructuredGridParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(IncToStructuredGridParameters, ApplicationParameters)
  typedef Increment<MODEL>               Increment_;
  typedef Geometry<MODEL>                Geometry_;

 public:
  typedef typename Increment_::ReadParameters_      IncrementParameters_;
  typedef typename Geometry_::Parameters_           GeometryParameters_;

  RequiredParameter<GeometryParameters_>        incGeometry{"increment geometry", this};
  RequiredParameter<Variables>                  vars{"variables", this};
  RequiredParameter<util::DateTime>             date{"date", this};
  RequiredParameter<IncrementParameters_>       increment{"increment", this};
  RequiredParameter<eckit::LocalConfiguration>  structuredGridInterp
                   {"structured grid interpolation", this};
};

template <typename MODEL> class StateEnsToStructuredGridParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(StateEnsToStructuredGridParameters, ApplicationParameters)
  typedef StateEnsembleParameters<MODEL> StateEnsembleParameters_;
  typedef Geometry<MODEL>                Geometry_;

 public:
  typedef typename Geometry_::Parameters_ GeometryParameters_;

  RequiredParameter<GeometryParameters_>       stateGeometry{"state geometry", this};
  RequiredParameter<StateEnsembleParameters_>  stateEnsemble{"states", this};
  RequiredParameter<eckit::LocalConfiguration> structuredGridInterp
                   {"structured grid interpolation", this};
};

template <typename MODEL> class IncEnsToStructuredGridParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(IncEnsToStructuredGridParameters, ApplicationParameters)
  typedef IncrementEnsembleParameters<MODEL> IncrementEnsembleParameters_;
  typedef Geometry<MODEL>                Geometry_;

 public:
  typedef typename Geometry_::Parameters_ GeometryParameters_;

  RequiredParameter<GeometryParameters_>          incrementGeometry{"increment geometry", this};
  RequiredParameter<IncrementEnsembleParameters_> incrementEnsemble{"increments", this};
  RequiredParameter<Variables>                    incrementVariables{"increment variables", this};
  RequiredParameter<eckit::LocalConfiguration>    structuredGridInterp
                   {"structured grid interpolation", this};
};

template <typename MODEL> class ConvertToStructuredGridParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(ConvertToStructuredGridParameters, ApplicationParameters)
  typedef StateEnsToStructuredGridParameters<MODEL>       StateEnsToStructuredGridParameters_;
  typedef IncEnsToStructuredGridParameters<MODEL>         IncEnsToStructuredGridParameters_;
  typedef StateToStructuredGridParameters<MODEL>          StateToStructuredGridParameters_;
  typedef IncToStructuredGridParameters<MODEL>            IncToStructuredGridParameters_;

 public:
  OptionalParameter<StateEnsToStructuredGridParameters_>
                   stateEnsToStructuredGrid{"state ensemble to structured grid", this};
  OptionalParameter<IncEnsToStructuredGridParameters_>
                   incEnsToStructuredGrid{"increment ensemble to structured grid", this};
  OptionalParameter<std::vector<StateToStructuredGridParameters_>>
                   stateToStructuredGrid{"states to structured grid", this};
  OptionalParameter<std::vector<IncToStructuredGridParameters_>>
                   incToStructuredGrid{"increments to structured grid", this};
};

template <typename MODEL> class ConvertToStructuredGrid : public Application {
  typedef Geometry<MODEL>                           Geometry_;
  typedef State<MODEL>                              State_;
  typedef StateEnsemble<MODEL>                      StateEnsemble_;
  typedef Increment<MODEL>                          Increment_;
  typedef IncrementEnsemble<MODEL>                  IncrementEnsemble_;
  typedef std::shared_ptr<IncrementEnsemble<MODEL>> EnsemblePtr_;
  typedef StructuredGridWriter<MODEL>                   StructuredGridGridWriter_;

  typedef ConvertToStructuredGridParameters<MODEL>          ConvertToStructuredGridParameters_;

 public:
// -----------------------------------------------------------------------------
  explicit ConvertToStructuredGrid(const eckit::mpi::Comm & comm = oops::mpi::world()) :
                                                              Application(comm) {}
// -----------------------------------------------------------------------------
  virtual ~ConvertToStructuredGrid() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
//  Deserialize parameters
    ConvertToStructuredGridParameters_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

// -----------------------------------------------------------------------------

//  Interpolate state ensemble if provided
    if (params.stateEnsToStructuredGrid.value() != boost::none) {
      Log::info() << "Interpolating State Ensemble" << std::endl;
      Geometry_ resol_(params.stateEnsToStructuredGrid.value()->stateGeometry, this->getComm());
      StateEnsemble_ statesToInterp_(resol_, params.stateEnsToStructuredGrid.value()->
                                                                 stateEnsemble.value());
      const eckit::LocalConfiguration structuredgridConf =
        params.stateEnsToStructuredGrid.value()->structuredGridInterp.value();
      const StructuredGridGridWriter_ structuredGridWriter_(structuredgridConf, resol_);
      size_t numstates = statesToInterp_.size();
      for (size_t jm=0; jm < numstates; jm++) {
        structuredGridWriter_.interpolateAndWrite(statesToInterp_[jm]);
        Log::test() << structuredGridWriter_ << std::endl;
      }
    }

// -----------------------------------------------------------------------------

//  Interpolate individual states if provided
    if (params.stateToStructuredGrid.value() != boost::none) {
      Log::info() << "Interpolating Individual State(s) " << std::endl;
      size_t numstates = params.stateToStructuredGrid.value()->size();
      for (size_t jm=0; jm < numstates; jm++) {
        Geometry_ resol_(params.stateToStructuredGrid.value()->at(jm).stateGeometry,
                         this->getComm());
        State_ stateToInterp_(resol_, params.stateToStructuredGrid.value()->at(jm).state.value());
        const eckit::LocalConfiguration structuredgridConf =
          params.stateToStructuredGrid.value()->at(jm).structuredGridInterp.value();
        const StructuredGridGridWriter_ structuredGridWriter_(structuredgridConf, resol_);
        structuredGridWriter_.interpolateAndWrite(stateToInterp_);
        Log::test() << structuredGridWriter_ << std::endl;
      }
    }

// -----------------------------------------------------------------------------

//  Interpolate individual increments if provided
    if (params.incToStructuredGrid.value() != boost::none) {
      Log::info() << "Interpolating Individual Increment(s)" << std::endl;
      size_t numstates = params.incToStructuredGrid.value()->size();
      for (size_t jm=0; jm < numstates; jm++) {
        Geometry_ resol_(params.incToStructuredGrid.value()->at(jm).incGeometry, this->getComm());
        Increment_ incToInterp_(resol_, params.incToStructuredGrid.value()->at(jm).vars.value(),
                                        params.incToStructuredGrid.value()->at(jm).date.value());
        incToInterp_.read(params.incToStructuredGrid.value()->at(jm).increment.value());
        const eckit::LocalConfiguration structuredgridConf =
          params.stateToStructuredGrid.value()->at(jm).structuredGridInterp.value();
        const StructuredGridGridWriter_ structuredGridWriter_(structuredgridConf, resol_);
        // This supports output on model levels only; to output on pressure levels would need to
        // read in a reference background from which to read the vertical pressure coordinate
        structuredGridWriter_.interpolateAndWrite(incToInterp_);
        Log::test() << structuredGridWriter_ << std::endl;
      }
    }
// -----------------------------------------------------------------------------

    if (params.incEnsToStructuredGrid.value() != boost::none) {
      Log::info() << "Interpolating Increment Ensemble" << std::endl;
      Geometry_ resol_(params.incEnsToStructuredGrid.value()->incrementGeometry, this->getComm());
      IncrementEnsemble_ incrementsToInterp_(resol_,
                            params.incEnsToStructuredGrid.value()->incrementVariables.value(),
                            params.incEnsToStructuredGrid.value()->incrementEnsemble.value());
      const eckit::LocalConfiguration structuredgridConf =
        params.incEnsToStructuredGrid.value()->structuredGridInterp.value();
      const StructuredGridGridWriter_ structuredGridWriter_(structuredgridConf, resol_);
      size_t numstates = incrementsToInterp_.size();
      for (size_t jm=0; jm < numstates; jm++) {
        // This supports output on model levels only; to output on pressure levels would need to
        // read in a reference background from which to read the vertical pressure coordinate
        structuredGridWriter_.interpolateAndWrite(incrementsToInterp_[jm]);
        Log::test() << structuredGridWriter_ << std::endl;
      }
    }

    return 0;
  }
// -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    ConvertToStructuredGridParameters_ params;
    params.outputSchema(outputPath);
  }
// -----------------------------------------------------------------------------
  void validateConfig(const eckit::Configuration & fullConfig) const override {
    ConvertToStructuredGridParameters_ params;
    params.validate(fullConfig);
  }
// -------------------------------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::ConvertToStructuredGrid<" + MODEL::name() + ">";
  }
};

}  // namespace oops

#endif  // OOPS_RUNS_CONVERTTOSTRUCTUREDGRID_H_
