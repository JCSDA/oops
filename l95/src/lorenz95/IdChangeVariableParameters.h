/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_IDCHANGEVARIABLEPARAMETERS_H_
#define LORENZ95_IDCHANGEVARIABLEPARAMETERS_H_

#include <ostream>
#include <string>

#include "oops/base/VariableChangeParametersBase.h"
#include "oops/util/Printable.h"

namespace lorenz95 {

// -------------------------------------------------------------------------------------------------
/// No change of variable parameters

class IdChangeVariableParameters : public oops::VariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(IdChangeVariableParameters, VariableChangeParametersBase)
 public:
  // No variable change. No additional parameters
};

// -------------------------------------------------------------------------------------------------

}  // namespace lorenz95
#endif  // LORENZ95_IDCHANGEVARIABLEPARAMETERS_H_
