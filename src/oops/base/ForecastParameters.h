/*
 * (C) Copyright 2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_FORECASTPARAMETERS_H_
#define OOPS_BASE_FORECASTPARAMETERS_H_

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Geometry.h"
#include "oops/base/LatLonGridPostProcessor.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Forecast parameters
template <typename MODEL> class ForecastParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(ForecastParameters, Parameters)

 public:
  typedef typename Geometry<MODEL>::Parameters_    GeometryParameters_;

  /// Geometry parameters.
  RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Model parameters.
  RequiredParameter<eckit::LocalConfiguration> model{"model", this};

  /// Initial state parameters.
  RequiredParameter<eckit::LocalConfiguration> initialCondition{"initial condition", this};

  /// Augmented model state.
  Parameter<eckit::LocalConfiguration> modelAuxControl{
    "model aux control", eckit::LocalConfiguration(), this};

  /// Forecast length.
  RequiredParameter<util::Duration> forecastLength{"forecast length", this};

  /// Where to write the output.
  RequiredParameter<eckit::LocalConfiguration> output{"output", this};
  OptionalParameter<LatLonGridPostProcessorParameters> latlonGridOutput{"forecast to latlon", this};

  /// Options passed to the object writing out forecast fields.
  Parameter<eckit::LocalConfiguration> prints{"prints", eckit::LocalConfiguration(), this};
};

// -----------------------------------------------------------------------------

/// Forecast TL parameters
template <typename MODEL> class ForecastTLParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(ForecastTLParameters, Parameters)

 public:
  typedef typename Geometry<MODEL>::Parameters_ GeometryParameters_;

  /// Geometry parameters.
  RequiredParameter<GeometryParameters_> geometry{"geometry", this};

  /// Linear model parameters.
  RequiredParameter<eckit::LocalConfiguration> model{"linear model", this};
};

// -----------------------------------------------------------------------------

/// Forecast AD parameters
template <typename MODEL> class ForecastADParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(ForecastADParameters, Parameters)

 public:
  typedef typename Geometry<MODEL>::Parameters_ GeometryParameters_;

  /// Adjoint initial condition output parameters.
  OptionalParameter<eckit::LocalConfiguration>
    adjointInitialConditionOutput{"adjoint initial condition output", this};

  /// Adjoint forecast output parameters.
  RequiredParameter<eckit::LocalConfiguration>
    adjointForecastOutput{"adjoint forecast output", this};

  /// Augmented model increment parameters.
  Parameter<eckit::LocalConfiguration> modelAuxIncrement{
    "model aux increment", eckit::LocalConfiguration(), this};
};

// -----------------------------------------------------------------------------

/// ForecastAspect parameters
template <typename MODEL> class ForecastAspectParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(ForecastAspectParameters, Parameters)

 public:
  typedef typename Geometry<MODEL>::Parameters_ GeometryParameters_;

  /// Verification resolution parameters.
  RequiredParameter<GeometryParameters_> verificationResolution{"verification resolution", this};

  /// Verification state parameters.
  RequiredParameter<eckit::LocalConfiguration> verificationConfig{"verification state", this};

  /// Adjoint initial condition norm parameters.
  RequiredParameter<eckit::LocalConfiguration> normConfig{"norm", this};
};

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_BASE_FORECASTPARAMETERS_H_
