/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_VARIABLECHANGEPARAMETERSBASE_H_
#define OOPS_BASE_VARIABLECHANGEPARAMETERSBASE_H_

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"

namespace oops {

/// \brief Base class of classes storing parameters controlling specific variable changes.
class VariableChangeParametersBase : public Parameters {
  OOPS_ABSTRACT_PARAMETERS(VariableChangeParametersBase, Parameters)
 public:
  OptionalParameter<Variables> inputVariables{"input variables", this};
  OptionalParameter<Variables> outputVariables{"output variables", this};
};

}  // namespace oops

#endif  // OOPS_BASE_VARIABLECHANGEPARAMETERSBASE_H_
