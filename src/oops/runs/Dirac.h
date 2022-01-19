/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_RUNS_DIRAC_H_
#define OOPS_RUNS_DIRAC_H_

#include <memory>
#include <sstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/instantiateCovarFactory.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
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

/// \brief Top-level options taken by the Dirac application.
template <typename MODEL> class DiracParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(DiracParameters, ApplicationParameters)

 public:
  typedef ModelSpaceCovarianceParametersWrapper<MODEL> CovarianceParameters_;
  typedef typename Geometry<MODEL>::Parameters_        GeometryParameters_;
  typedef typename State<MODEL>::Parameters_           StateParameters_;
  typedef StateParametersND<MODEL>                     StateParametersND_;

  /// Geometry parameters.
  RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Initial state parameters.
  RequiredParameter<StateParametersND_> initialCondition{"initial condition", this};

  /// Background error covariance model.
  RequiredParameter<CovarianceParameters_> backgroundError{"background error", this};

  /// Dirac location/variables parameters.
  RequiredParameter<eckit::LocalConfiguration> dirac{"dirac", this};

  /// Where to write the output(s) of Dirac tests
  RequiredParameter<eckit::LocalConfiguration> outputDirac{"output dirac", this};

  /// Where to write the output of randomized variance.
  OptionalParameter<eckit::LocalConfiguration> outputVariance{"output variance", this};
};

// -----------------------------------------------------------------------------

template <typename MODEL> class Dirac : public Application {
  typedef ModelSpaceCovarianceBase<MODEL>           CovarianceBase_;
  typedef CovarianceFactory<MODEL>                  CovarianceFactory_;
  typedef ModelSpaceCovarianceParametersBase<MODEL> CovarianceParametersBase_;
  typedef Geometry<MODEL>                           Geometry_;
  typedef Increment<MODEL>                          Increment_;
  typedef State<MODEL>                              State_;
  typedef Localization<MODEL>                       Localization_;
  typedef IncrementEnsemble<MODEL>                  Ensemble_;
  typedef std::shared_ptr<IncrementEnsemble<MODEL>> EnsemblePtr_;

  typedef DiracParameters<MODEL>                    DiracParameters_;

 public:
// -----------------------------------------------------------------------------
  explicit Dirac(const eckit::mpi::Comm & comm = oops::mpi::world()) : Application(comm) {
    instantiateCovarFactory<MODEL>();
  }
// -----------------------------------------------------------------------------
  virtual ~Dirac() {}
// -----------------------------------------------------------------------------
  int execute(const eckit::Configuration & fullConfig) const override {
//  Deserialize parameters
    DiracParameters_ params;
    params.validateAndDeserialize(fullConfig);

//  Define number of subwindows
    const eckit::LocalConfiguration backgroundConfig(fullConfig, "initial condition");
    std::vector<eckit::LocalConfiguration> confs;
    backgroundConfig.get("states", confs);
    size_t nslots = confs.size();

//  Define space and time communicators
    const eckit::mpi::Comm * commSpace = &this->getComm();
    const eckit::mpi::Comm * commTime = &oops::mpi::myself();
    if (nslots > 1) {
      size_t ntasks = this->getComm().size();
      ASSERT(ntasks % nslots == 0);
      size_t myrank = this->getComm().rank();
      size_t ntaskpslot = ntasks / nslots;
      size_t myslot = myrank / ntaskpslot;

      // Create a communicator for same sub-window, to be used for communications in space
      std::string sgeom = "comm_geom_" + std::to_string(myslot);
      char const *geomName = sgeom.c_str();
      commSpace = &this->getComm().split(myslot, geomName);
      ASSERT(commSpace->size() == ntaskpslot);

      // Create a communicator for same local area, to be used for communications in time
      size_t myarea = commSpace->rank();
      std::string stime = "comm_time_" + std::to_string(myarea);
      char const *timeName = stime.c_str();
      commTime = &this->getComm().split(myarea, timeName);
      ASSERT(commTime->size() == nslots);
    }

//  Setup resolution
    const Geometry_ resol(params.geometry, *commSpace, *commTime);

//  Setup initial state
    const State_ xx(resol, params.initialCondition);

//  Setup variables
    const Variables vars = xx.variables();

//  Setup time
    util::DateTime time = xx.validTime();

//  Setup Dirac field
    Increment_ dxi(resol, vars, time);
    dxi.dirac(params.dirac);
    Log::test() << "Input Dirac increment: " << dxi << std::endl;

//  Go recursively through the covariance configuration
    const CovarianceParametersBase_ &covarParams =
        params.backgroundError.value().covarianceParameters;
    eckit::LocalConfiguration covarConf(covarParams.toConfiguration());
    std::string id;
    dirac(covarConf, params.outputDirac.value(), id, resol, xx, dxi);

//  Variance randomization
    const boost::optional<eckit::LocalConfiguration> &outputVariance =
        params.outputVariance.value();
    if (outputVariance != boost::none) {
      // Setup variance
      Increment_ dx(resol, vars, time);
      Increment_ dxsq(resol, vars, time);
      Increment_ mean(resol, vars, time);
      Increment_ variance(resol, vars, time);
      mean.zero();
      variance.zero();

      // Covariance
      std::unique_ptr<CovarianceBase_> Bmat(CovarianceFactory_::create(
          covarConf, resol, vars, xx, xx));

      // Randomization
      for (size_t ie = 0; ie < Bmat->randomizationSize(); ++ie) {
        Bmat->randomize(dx);
        dx -= mean;
        dxsq = dx;
        dxsq.schur_product_with(dx);
        double rk_var = static_cast<double>(ie)/static_cast<double>(ie+1);
        double rk_mean = 1.0/static_cast<double>(ie+1);
        variance.axpy(rk_var, dxsq, false);
        mean.axpy(rk_mean, dx, false);
      }
      double rk_norm = 1.0/static_cast<double>(Bmat->randomizationSize()-1);
      variance *= rk_norm;

      // Write increment
      variance.write(*(params.outputVariance.value()));
      Log::test() << "Randomized variance: " << variance << std::endl;
    }

    return 0;
  }
// -----------------------------------------------------------------------------
  void outputSchema(const std::string & outputPath) const override {
    DiracParameters_ params;
    params.outputSchema(outputPath);
  }
// -----------------------------------------------------------------------------
 private:
  std::string appname() const override {
    return "oops::Dirac<" + MODEL::name() + ">";
  }

  void dirac(const eckit::LocalConfiguration & covarConfig,
             const eckit::LocalConfiguration & outputConfig,
             std::string & id,
             const Geometry_ & resol, const State_ & xx,
             const Increment_ & dxi) const {
    Log::debug() << "Input ID is " << id << std::endl;
    // Define output increment
    Increment_ dxo(dxi, false);

    // Covariance
    std::unique_ptr<CovarianceBase_> Bmat(CovarianceFactory_::create(
        covarConfig, resol, xx.variables(), xx, xx));

    // Multiply
    Bmat->multiply(dxi, dxo);

    // Copy configuration
    eckit::LocalConfiguration outputBConf(outputConfig);

    // Update ID
    if (id != "") id.append("_");
    id.append(Bmat->covarianceModel());

    // Seek and replace %id% with id, recursively
    util::seekAndReplace(outputBConf, "%id%", id);

    // Write output increment
    dxo.write(outputBConf);
    Log::test() << "Covariance(" << id << ") * Increment: " << dxo << std::endl;

    // Look for hybrid or ensemble covariance models
    const std::string covarianceModel(covarConfig.getString("covariance model"));
    if (covarianceModel == "hybrid") {
      std::vector<eckit::LocalConfiguration> confs;
      covarConfig.get("components", confs);
      for (const auto & conf : confs) {
        std::string idC(id);
        const eckit::LocalConfiguration componentConfig(conf, "covariance");
        dirac(componentConfig, outputConfig, idC, resol, xx, dxi);
      }
    }
    if (covarianceModel == "ensemble" && covarConfig.has("localization")) {
      // Localization configuration
      eckit::LocalConfiguration locConfig(covarConfig.getSubConfiguration("localization"));

      // Define output increment
      Increment_ dxo(dxi);

      // Setup localization
      Localization_ Lmat(resol, locConfig);

      // Apply localization
      Lmat.multiply(dxo);

      // Copy configuration
      eckit::LocalConfiguration outputLConf(outputConfig);

      // Update ID
      std::string idL(id);
      idL.append("_localization");

      // Seek and replace %id% with id, recursively
      util::seekAndReplace(outputLConf, "%id%", idL);

      // Write output increment
      dxo.write(outputLConf);
      Log::test() << "Localization(" << id << ") * Increment: " << dxo << std::endl;
    }
  }
// -----------------------------------------------------------------------------
};

}  // namespace oops
#endif  // OOPS_RUNS_DIRAC_H_
