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
  typedef typename State<MODEL>::WriteParameters_         StateWriteParameters_;
  typedef StateEnsembleParameters<MODEL>                  StateEnsembleParameters_;

  /// Geometry parameters.
  RequiredParameter<GeometryParameters_> resolConfig{"geometry", this};

  /// Ensemble parameters.
  RequiredParameter<StateEnsembleParameters_> ensembleConfig{"ensemble", this};

  /// Output increment parameters for variance.
  OptionalParameter<IncrementWriteParameters_> outputStdDevConfig{"standard deviation output",
                                                this};
  OptionalParameter<LatLonGridWriterParameters> outputStdDevConfigLL{"standard deviation to latlon",
                                                 this};
  OptionalParameter<IncrementWriteParameters_> outputVarConfig{"variance output", this};
  OptionalParameter<LatLonGridWriterParameters> outputVarConfigLL{"ensvariance to latlon",
                                                 this};

  /// Output state parameters for mean.
  OptionalParameter<StateWriteParameters_> outputMeanConfig{"mean output", this};
  OptionalParameter<LatLonGridWriterParameters> outputMeanConfigLL{"ensmean to latlon",
                                                 this};
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
    if (params.outputMeanConfig.value() != boost::none)
        ensmean.write(params.outputMeanConfig.value().value());

    if (params.outputMeanConfigLL.value() != boost::none) {
      const LatLonGridWriter<MODEL> latlon(params.outputMeanConfigLL.value().value(), resol);
      latlon.interpolateAndWrite(ensmean);
    }
    Log::test() << "Mean: " << std::endl << ensmean << std::endl;

//  Write variance to file
    if (params.outputVarConfig.value() != boost::none)
        sigb2.write(params.outputVarConfig.value().value());

    if (params.outputVarConfigLL.value() != boost::none) {
      const LatLonGridWriter<MODEL> latlon(params.outputVarConfigLL.value().value(), resol);
      latlon.interpolateAndWrite(sigb2, ensmean);
    }
    Log::test() << "Variance: " << std::endl << sigb2 << std::endl;

//  Write standard deviation to file
    if (params.outputStdDevConfig.value() != boost::none)
        sigb.write(params.outputStdDevConfig.value().value());

    if (params.outputStdDevConfigLL.value() != boost::none) {
      const LatLonGridWriter<MODEL> latlon(params.outputStdDevConfigLL.value().value(), resol);
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
