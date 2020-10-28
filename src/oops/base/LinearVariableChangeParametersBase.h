/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_LINEARVARIABLECHANGEPARAMETERSBASE_H_
#define OOPS_BASE_LINEARVARIABLECHANGEPARAMETERSBASE_H_

#include <string>

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/Variables.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"

namespace oops {

/// \brief Base class of classes storing parameters controlling specific linear variable changes.
class LinearVariableChangeParametersBase : public Parameters {
  OOPS_ABSTRACT_PARAMETERS(LinearVariableChangeParametersBase, Parameters)
 public:
  /// \brief Variable change type.
  ///
  /// \note This parameter is marked as optional because it is only required in certain
  /// circumstances (e.g. when variable change parameters are deserialized into a
  /// LinearVariableChangeParametersWrapper and used by LinearVariableChangeFactory to instantiate
  /// a tangent linear variable change whose type is determined at runtime), but not others (e.g.
  /// in tests written with a particular variable change in mind).
  /// LinearVariableChangeParametersWrapper will throw an exception if this parameter is not
  /// provided.
  OptionalParameter<std::string> variableChange{"variable change", this};

  OptionalParameter<Variables> inputVariables{"input variables", this};
  OptionalParameter<Variables> outputVariables{"output variables", this};
};

}  // namespace oops

#endif  // OOPS_BASE_LINEARVARIABLECHANGEPARAMETERSBASE_H_
