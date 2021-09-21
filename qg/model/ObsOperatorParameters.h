/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QG_MODEL_OBSOPERATORPARAMETERS_H_
#define QG_MODEL_OBSOPERATORPARAMETERS_H_

#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/Parameters.h"

namespace qg {

/// Parameters controlling the observation operator for the QG model.
class ObservationParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObservationParameters, Parameters);
 public:
  oops::ConfigurationParameter config{this};
};

// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_OBSOPERATORPARAMETERS_H_
