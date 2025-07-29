/*
 * (C) Crown copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_ENSEMBLEINFLATION_H_
#define OOPS_RUNS_ENSEMBLEINFLATION_H_

#include <memory>
#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/IncrementSet.h"
#include "oops/base/Inflation.h"
#include "oops/base/InflationBase.h"
#include "oops/base/instantiateInflationFactory.h"
#include "oops/base/State.h"
#include "oops/base/StateSet.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"

namespace oops {

template <typename MODEL> class EnsembleInflation : public Application {
  typedef Geometry<MODEL>                   Geometry_;
  typedef Increment<MODEL>                  Increment_;
  typedef IncrementSet<MODEL>               IncrementSet_;
  typedef State<MODEL>                      State_;
  typedef StateSet<MODEL>                   StateSet_;

 public:
// -----------------------------------------------------------------------------

  explicit EnsembleInflation(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {
    instantiateInflationFactory<MODEL>();
  }

// -----------------------------------------------------------------------------

  virtual ~EnsembleInflation() = default;

// -----------------------------------------------------------------------------

  int execute(const eckit::Configuration & fullConfig) const override {
    // Setup geometry
    const Geometry_ geometry(eckit::LocalConfiguration(fullConfig, "geometry"),
                             this->getComm(), oops::mpi::myself());

    // Read all background ensemble members
    StateSet_ bgens(geometry, eckit::LocalConfiguration(fullConfig, "background"));

    // Get inflation subconfigurations
    eckit::LocalConfiguration infConf(fullConfig, "inflation");
    std::vector<eckit::LocalConfiguration> subconfigs = infConf.getSubConfigurations();

    // Carry out inflation depending on whether analysis is in the form of increments or states
    eckit::LocalConfiguration anConf(fullConfig, "analysis");
    if (fullConfig.getString("analysis type") == "state") {
      Inflation<MODEL, StateSet_> inflation(anConf, geometry, bgens);
      inflation.calculate(subconfigs);
      inflation.save(eckit::LocalConfiguration(fullConfig, "output"));
    } else {
      // Setup analysis variables
      Variables anvars(fullConfig, "analysis variables");
      Inflation<MODEL, IncrementSet_> inflation(anConf, geometry, bgens, anvars);
      inflation.calculate(subconfigs);
      inflation.save(eckit::LocalConfiguration(fullConfig, "output"));
    }
    return 0;
}
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::EnsembleInflation<" + MODEL::name() + ">";
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_ENSEMBLEINFLATION_H_
