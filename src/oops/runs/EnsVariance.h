/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_RUNS_ENSVARIANCE_H_
#define OOPS_RUNS_ENSVARIANCE_H_

#include <memory>
#include <string>
#include <vector>


#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Geometry.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the EnsVariance application.
template <typename MODEL>
class EnsVarianceParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(EnsVarianceParameters, ApplicationParameters)

 public:
  typedef typename Geometry<MODEL>::Parameters_ GeometryParameters_;
  typedef typename State<MODEL>::Parameters_ StateParameters_;
  typedef typename Increment<MODEL>::WriteParameters_ IncrementWriteParameters_;

  /// Geometry parameters.
  RequiredParameter<GeometryParameters_> resolConfig{"geometry", this};

  /// Background parameters.
  RequiredParameter<StateParameters_> bkgConfig{"background", this};

  /// Ensemble parameters.
  RequiredParameter<eckit::LocalConfiguration> ensembleConfig{"ensemble", this};

  /// Output increment parameters.
  RequiredParameter<IncrementWriteParameters_> outputConfig{"variance output", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class EnsVariance : public Application {
  typedef IncrementEnsemble<MODEL>                 Ensemble_;
  typedef Geometry<MODEL>                          Geometry_;
  typedef Increment<MODEL>                         Increment_;
  typedef State<MODEL>                             State_;
  typedef EnsVarianceParameters<MODEL>             EnsVarianceParameters_;

 public:
  // -----------------------------------------------------------------------------
  explicit EnsVariance(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {}
  // -----------------------------------------------------------------------------
  virtual ~EnsVariance() {}
  // -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
//  Deserialize parameters
    EnsVarianceParameters_ params;
    params.validateAndDeserialize(fullConfig);

//  Setup Geometry
    const Geometry_ resol(params.resolConfig, this->getComm());

//  Setup background
    State_ xx(resol, params.bkgConfig);

//  Compute transformed ensemble perturbations
//         ens_k = K^-1 dx_k
    const eckit::LocalConfiguration ensConfig = params.ensembleConfig;
    Variables vars(ensConfig, "output variables");
    Ensemble_ ens_k(ensConfig, xx, xx, resol, vars);

//  Get ensemble size
    unsigned nm = ens_k.size();

//  Compute ensemble standard deviation
    Increment_ km1dx(ens_k[0]);
    km1dx.zero();
    Increment_ sigb2(km1dx);
    sigb2.zero();

    for (unsigned jj = 0; jj < nm; ++jj) {
      km1dx = ens_k[jj];

//    Accumulate km1dx^2
      km1dx.schur_product_with(km1dx);
      sigb2 += km1dx;
    }
    const double rk = 1.0/(static_cast<double>(nm) - 1.0);
    sigb2 *= rk;

//  Write variance to file
    sigb2.write(params.outputConfig);
    Log::test() << "Variance: " << std::endl << sigb2 << std::endl;

    return 0;
  }
  // -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    EnsVarianceParameters_ params;
    params.outputSchema(outputPath);
  }
  // -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::EnsVariance<" + MODEL::name() + ">";
  }
  // -----------------------------------------------------------------------------
};

}  // namespace oops

#endif  // OOPS_RUNS_ENSVARIANCE_H_
