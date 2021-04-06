/*
 * (C) Crown copyright 2021, Met Office.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef LORENZ95_OBSBIASPARAMETERS_H_
#define LORENZ95_OBSBIASPARAMETERS_H_

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"

namespace lorenz95 {

/// Parameters taken by the ObsBias, ObsBiasCorrection and ObsBiasCovariance classes.

// -----------------------------------------------------------------------------

class ObsBiasCovarianceParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsBiasCovarianceParameters, Parameters)

 public:
  oops::OptionalParameter<double> standardDeviation{"standard_deviation", this};
};

// -----------------------------------------------------------------------------

class ObsBiasParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsBiasParameters, Parameters)

 public:
  oops::OptionalParameter<double> bias{"bias", this};
  oops::OptionalParameter<ObsBiasCovarianceParameters> covariance{"covariance", this};
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_OBSBIASPARAMETERS_H_
