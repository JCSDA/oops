/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_IDCHANGEVARTLADL95PARAMS_H_
#define LORENZ95_IDCHANGEVARTLADL95PARAMS_H_

#include <ostream>
#include <string>

#include "oops/base/LinearVariableChangeParametersBase.h"

namespace lorenz95 {

// -------------------------------------------------------------------------------------------------
/// No change of linear variable parameters

class IdChangeVarTLADL95Params : public oops::LinearVariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(IdChangeVarTLADL95Params, LinearVariableChangeParametersBase)
 public:
  // No linear variable change. No additional parameters
};

// -------------------------------------------------------------------------------------------------

}  // namespace lorenz95
#endif  // LORENZ95_IDCHANGEVARTLADL95PARAMS_H_
