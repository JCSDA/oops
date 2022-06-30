/*
 * (C) Copyright 2020 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_MODELBIASCOVARIANCEPARAMETERS_H_
#define LORENZ95_MODELBIASCOVARIANCEPARAMETERS_H_

#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"

namespace lorenz95 {

/// \brief Parameters controlling the bias covariance and increment.
class ModelBiasCovarianceParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ModelBiasCovarianceParameters, Parameters)

 public:
  oops::OptionalParameter<double> standardDeviation{"standard_deviation", this};
};

}  // namespace lorenz95

#endif  // LORENZ95_MODELBIASCOVARIANCEPARAMETERS_H_
