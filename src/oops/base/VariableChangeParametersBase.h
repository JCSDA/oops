/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_VARIABLECHANGEPARAMETERSBASE_H_
#define OOPS_BASE_VARIABLECHANGEPARAMETERSBASE_H_

#include <string>

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/ConfigurationParameter.h"
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

/// \brief A subclass of VariableChangeParametersBase storing the values of all options in a
/// single Configuration object.
///
/// This object can be accessed by calling the value() method of the \p config member variable.
///
/// The ConfigurationParameter class does not perform any parameter validation; models using
/// GenericVariableChangeParameters should therefore ideally be refactored, replacing this
/// class with a dedicated subclass of VariableChangeParametersBase storing each parameter in
/// a separate (Optional/Required)Parameter object.
class GenericVariableChangeParameters : public VariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(GenericVariableChangeParameters, VariableChangeParametersBase)
 public:
  ConfigurationParameter config{this};
};


}  // namespace oops

#endif  // OOPS_BASE_VARIABLECHANGEPARAMETERSBASE_H_
