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

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL> class ConvertToStructuredGrid : public Application {
  typedef Geometry<MODEL>                           Geometry_;
  typedef State<MODEL>                              State_;
  typedef StateEnsemble<MODEL>                      StateEnsemble_;
  typedef Increment<MODEL>                          Increment_;
  typedef IncrementEnsemble<MODEL>                  IncrementEnsemble_;
  typedef std::shared_ptr<IncrementEnsemble<MODEL>> EnsemblePtr_;
  typedef StructuredGridWriter<MODEL>               StructuredGridGridWriter_;

 public:
// -----------------------------------------------------------------------------
  explicit ConvertToStructuredGrid(const eckit::mpi::Comm & comm = oops::mpi::world()) :
                                                              Application(comm) {}
// -----------------------------------------------------------------------------
  virtual ~ConvertToStructuredGrid() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
//  Interpolate state ensemble if provided
    const eckit::LocalConfiguration stateEns =
        fullConfig.getSubConfiguration("state ensemble to structured grid");
    if (!stateEns.empty()) {
      Log::info() << "Interpolating State Ensemble" << std::endl;
      Geometry_ resol_(eckit::LocalConfiguration(stateEns, "state geometry"), this->getComm());
      eckit::LocalConfiguration statesConf(stateEns, "states");
      StateEnsemble_ statesToInterp_(resol_, statesConf);
      const eckit::LocalConfiguration structuredgridConf(stateEns, "structured grid interpolation");
      const StructuredGridGridWriter_ structuredGridWriter_(structuredgridConf, resol_);
      size_t numstates = statesToInterp_.size();
      for (size_t jm=0; jm < numstates; jm++) {
        structuredGridWriter_.interpolateAndWrite(statesToInterp_[jm]);
        Log::test() << structuredGridWriter_ << std::endl;
      }
    }

// -----------------------------------------------------------------------------

//  Interpolate individual states if provided
    if (fullConfig.has("states to structured grid")) {
      const std::vector<eckit::LocalConfiguration> statesConf =
          fullConfig.getSubConfigurations("states to structured grid");
      for (size_t jm=0; jm < statesConf.size(); jm++) {
        Geometry_ resol_(eckit::LocalConfiguration(statesConf.at(jm), "state geometry"),
                         this->getComm());
        State_ stateToInterp_(resol_, eckit::LocalConfiguration(statesConf.at(jm), "state"));
        const eckit::LocalConfiguration structuredgridConf(statesConf.at(jm),
                                                           "structured grid interpolation");
        const StructuredGridGridWriter_ structuredGridWriter_(structuredgridConf, resol_);
        structuredGridWriter_.interpolateAndWrite(stateToInterp_);
        Log::test() << structuredGridWriter_ << std::endl;
      }
    }

// -----------------------------------------------------------------------------

//  Interpolate individual increments if provided
    if (fullConfig.has("increments to structured grid")) {
      const std::vector<eckit::LocalConfiguration> incsConf =
          fullConfig.getSubConfigurations("increments to structured grid");
      for (size_t jm=0; jm < incsConf.size(); jm++) {
        Geometry_ resol_(eckit::LocalConfiguration(incsConf.at(jm), "increment geometry"),
                         this->getComm());
        const Variables vars(incsConf[jm], "variables");
        const util::DateTime tt(incsConf[jm].getString("date"));
        Increment_ incToInterp_(resol_, vars, tt);
        incToInterp_.read(eckit::LocalConfiguration(incsConf.at(jm), "increment"));
        const eckit::LocalConfiguration structuredgridConf(incsConf.at(jm),
                                                           "structured grid interpolation");
        const StructuredGridGridWriter_ structuredGridWriter_(structuredgridConf, resol_);
        // This supports output on model levels only; to output on pressure levels would need to
        // read in a reference background from which to read the vertical pressure coordinate
        structuredGridWriter_.interpolateAndWrite(incToInterp_);
        Log::test() << structuredGridWriter_ << std::endl;
      }
    }

// -----------------------------------------------------------------------------

    const eckit::LocalConfiguration incEns =
        fullConfig.getSubConfiguration("increment ensemble to structured grid");
    if (!incEns.empty()) {
      Log::info() << "Interpolating Increment Ensemble" << std::endl;
      Geometry_ resol_(eckit::LocalConfiguration(incEns, "increment geometry"), this->getComm());
      const Variables vars(incEns, "increment variables");
      eckit::LocalConfiguration incsConf(incEns, "increments");
      IncrementEnsemble_ incrementsToInterp_(resol_, vars, incsConf);
      const eckit::LocalConfiguration structuredgridConf(incEns, "structured grid interpolation");
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
 private:
  std::string appname() const override {
    return "oops::ConvertToStructuredGrid<" + MODEL::name() + ">";
  }
};

}  // namespace oops

#endif  // OOPS_RUNS_CONVERTTOSTRUCTUREDGRID_H_
