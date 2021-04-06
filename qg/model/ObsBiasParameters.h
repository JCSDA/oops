/*
 * (C) Crown copyright 2021, Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QG_MODEL_OBSBIASPARAMETERS_H_
#define QG_MODEL_OBSBIASPARAMETERS_H_

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"

namespace qg {

/// Parameters taken by the ObsBias, ObsBiasCorrection and ObsBiasCovariance classes.

// -----------------------------------------------------------------------------

class ObsBiasCovarianceParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsBiasCovarianceParameters, Parameters)

 public:
  oops::OptionalParameter<double> stream{"stream", this};
  oops::OptionalParameter<double> uwind{"uwind", this};
  oops::OptionalParameter<double> vwind{"vwind", this};
  oops::OptionalParameter<double> wspeed{"wspeed", this};
};

// -----------------------------------------------------------------------------

class ObsBiasParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsBiasParameters, Parameters)

 public:
  oops::OptionalParameter<double> stream{"stream", this};
  oops::OptionalParameter<double> uwind{"uwind", this};
  oops::OptionalParameter<double> vwind{"vwind", this};
  oops::OptionalParameter<double> wspeed{"wspeed", this};
  oops::OptionalParameter<ObsBiasCovarianceParameters> covariance{"covariance", this};
};

// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_OBSBIASPARAMETERS_H_
