/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QG_MODEL_CHANGEVARQGPARAMETERS_H_
#define QG_MODEL_CHANGEVARQGPARAMETERS_H_

#include <ostream>
#include <string>

#include "oops/base/VariableChangeParametersBase.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/Printable.h"

namespace qg {

// -------------------------------------------------------------------------------------------------
/// No change of variable parameters

class ChangeVarQGParameters : public oops::VariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(ChangeVarQGParameters, VariableChangeParametersBase)
 public:
};

// -------------------------------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_CHANGEVARQGPARAMETERS_H_
