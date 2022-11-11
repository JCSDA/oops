/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_ENSMEANANDVARIANCE_H_
#define OOPS_RUNS_ENSMEANANDVARIANCE_H_

#include <memory>
#include <string>
#include <vector>


#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/State.h"
#include "oops/base/StateEnsemble.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the EnsMeanAndVariance application.
template <typename MODEL>
class EnsMeanAndVarianceParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(EnsMeanAndVarianceParameters, ApplicationParameters)

 public:
  typedef typename Geometry<MODEL>::Parameters_        GeometryParameters_;
  typedef typename Increment<MODEL>::WriteParameters_  IncrementWriteParameters_;
  typedef typename State<MODEL>::WriteParameters_      StateWriteParameters_;
  typedef StateEnsembleParameters<MODEL>               StateEnsembleParameters_;

  /// Geometry parameters.
  RequiredParameter<GeometryParameters_> resolConfig{"geometry", this};

  /// Ensemble parameters.
  RequiredParameter<StateEnsembleParameters_> ensembleConfig{"ensemble", this};

  /// Output increment parameters for variance.
  RequiredParameter<IncrementWriteParameters_> outputVarConfig{"variance output", this};

  /// Output state parameters for mean.
  RequiredParameter<StateWriteParameters_> outputMeanConfig{"mean output", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class EnsMeanAndVariance : public Application {
  typedef Geometry<MODEL>                          Geometry_;
  typedef Increment<MODEL>                         Increment_;
  typedef State<MODEL>                             State_;
  typedef StateEnsemble<MODEL>                     StateEnsemble_;
  typedef EnsMeanAndVarianceParameters<MODEL>      EnsMeanAndVarianceParameters_;

 public:
  // -----------------------------------------------------------------------------
  explicit EnsMeanAndVariance(const eckit::mpi::Comm & comm = oops::mpi::world())
    : Application(comm) {}
  // -----------------------------------------------------------------------------
  virtual ~EnsMeanAndVariance() {}
  // -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig, bool validate) const override {
//  Deserialize parameters
    EnsMeanAndVarianceParameters_ params;
    if (validate) params.validate(fullConfig);
    params.deserialize(fullConfig);

//  Setup Geometry
    const Geometry_ resol(params.resolConfig, this->getComm());

//  Setup ensemble of states
    const StateEnsemble_ stateEnsemble(resol, params.ensembleConfig);
    const size_t nm = stateEnsemble.size();
    const State_ ensmean = stateEnsemble.mean();

//  Compute ensemble standard deviation
    Increment_ km1dx(resol, ensmean.variables(), ensmean.validTime());
    km1dx.zero();
    Increment_ sigb2(km1dx);
    sigb2.zero();

    for (unsigned jj = 0; jj < nm; ++jj) {
      km1dx.diff(stateEnsemble[jj], ensmean);

//    Accumulate km1dx^2
      km1dx.schur_product_with(km1dx);
      sigb2 += km1dx;
    }

    const double rk = 1.0/(static_cast<double>(nm) - 1.0);
    sigb2 *= rk;

//  Write mean to file
    ensmean.write(params.outputMeanConfig);
    Log::test() << "Mean: " << std::endl << ensmean << std::endl;

//  Write variance to file
    sigb2.write(params.outputVarConfig);
    Log::test() << "Variance: " << std::endl << sigb2 << std::endl;

    return 0;
  }
  // -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    EnsMeanAndVarianceParameters_ params;
    params.outputSchema(outputPath);
  }
// -----------------------------------------------------------------------------
  void validateConfig(const eckit::Configuration & fullConfig) const override {
    EnsMeanAndVarianceParameters_ params;
    params.validate(fullConfig);
  }
  // -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::EnsMeanAndVariance<" + MODEL::name() + ">";
  }
  // -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_ENSMEANANDVARIANCE_H_
