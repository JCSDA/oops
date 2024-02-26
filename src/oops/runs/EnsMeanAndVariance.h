/*
 * (C) Copyright 2022 UCAR
 * (C) Crown copyright 2024, Met Office
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
#include "oops/base/LatLonGridWriter.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/base/StateEnsemble.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/runs/Application.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Top-level options taken by the EnsMeanAndVariance application.
template <typename MODEL>
class EnsMeanAndVarianceParameters : public ApplicationParameters {
  OOPS_CONCRETE_PARAMETERS(EnsMeanAndVarianceParameters, ApplicationParameters)

 public:
  typedef typename Geometry<MODEL>::Parameters_           GeometryParameters_;
  typedef typename Increment<MODEL>::WriteParameters_     IncrementWriteParameters_;
  typedef StateEnsembleParameters<MODEL>                  StateEnsembleParameters_;

  /// Geometry parameters.
  RequiredParameter<GeometryParameters_> resolConfig{"geometry", this};

  /// Ensemble parameters.
  RequiredParameter<StateEnsembleParameters_> ensembleConfig{"ensemble", this};

  /// Output increment parameters for variance.
  OptionalParameter<IncrementWriteParameters_> outputStdDevConfig{"standard deviation output",
                                                 this};
  OptionalParameter<eckit::LocalConfiguration> outputStdDevConfigLL{"standard deviation to latlon",
                                                 this};
  OptionalParameter<IncrementWriteParameters_> outputVarConfig{"variance output", this};
  OptionalParameter<eckit::LocalConfiguration> outputVarConfigLL{"ensvariance to latlon", this};

  /// Output state parameters for mean.
  OptionalParameter<eckit::LocalConfiguration> outputMeanConfig{"mean output", this};
  OptionalParameter<eckit::LocalConfiguration> outputMeanConfigLL{"ensmean to latlon", this};
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
    const State_ ensmean = stateEnsemble.mean();
    const Increment_ sigb2 = stateEnsemble.variance();
    const Increment_ sigb = stateEnsemble.stddev();

//  Write mean to file
    if (fullConfig.has("mean output"))
      ensmean.write(eckit::LocalConfiguration(fullConfig, "mean output"));

    if (params.outputMeanConfigLL.value() != boost::none) {
      const eckit::LocalConfiguration latlonConf = params.outputMeanConfigLL.value().value();
      const LatLonGridWriter<MODEL> latlon(latlonConf, resol);
      latlon.interpolateAndWrite(ensmean);
    }
    Log::test() << "Mean: " << std::endl << ensmean << std::endl;

//  Write variance to file
    if (params.outputVarConfig.value() != boost::none)
        sigb2.write(params.outputVarConfig.value().value());

    if (params.outputVarConfigLL.value() != boost::none) {
      const eckit::LocalConfiguration latlonConf = params.outputVarConfigLL.value().value();
      const LatLonGridWriter<MODEL> latlon(latlonConf, resol);
      latlon.interpolateAndWrite(sigb2, ensmean);
    }
    Log::test() << "Variance: " << std::endl << sigb2 << std::endl;

//  Write standard deviation to file
    if (params.outputStdDevConfig.value() != boost::none)
        sigb.write(params.outputStdDevConfig.value().value());

    if (params.outputStdDevConfigLL.value() != boost::none) {
      const eckit::LocalConfiguration latlonConf = params.outputStdDevConfigLL.value().value();
      const LatLonGridWriter<MODEL> latlon(latlonConf, resol);
      latlon.interpolateAndWrite(sigb, ensmean);
    }
    Log::test() << "Standard Deviation: " << std::endl << sigb << std::endl;


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
