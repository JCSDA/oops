/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QG_MODEL_ERRORCOVARIANCEPARAMETERS_H_
#define QG_MODEL_ERRORCOVARIANCEPARAMETERS_H_

#include "oops/base/ModelSpaceCovarianceParametersBase.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace qg {
  class QgTraits;
// -----------------------------------------------------------------------------
/// \brief Parameters passed to the ErrorCovarianceParameters class.

class ErrorCovarianceParameters : public oops::ModelSpaceCovarianceParametersBase<qg::QgTraits>{
  OOPS_CONCRETE_PARAMETERS(ErrorCovarianceParameters,
                           oops::ModelSpaceCovarianceParametersBase<qg::QgTraits>)

 public:
  oops::RequiredParameter<double> standard_deviation{"standard_deviation", this};
  oops::RequiredParameter<double> horizontal_length_scale{"horizontal_length_scale", this};
  oops::RequiredParameter<double> vertical_length_scale{"vertical_length_scale", this};
  oops::RequiredParameter<double> maximum_condition_number{"maximum_condition_number", this};
  oops::Parameter<int> randomization_seed{"randomization_seed", 7, this};
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_ERRORCOVARIANCEPARAMETERS_H_
