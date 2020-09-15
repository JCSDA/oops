/*
 * (C) Copyright 2020 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_PARAMETERCONSTRAINT_H_
#define OOPS_UTIL_PARAMETERS_PARAMETERCONSTRAINT_H_

#include <memory>
#include <string>

#include "oops/util/parameters/ObjectJsonSchema.h"

namespace oops {

/// \brief Interface for classes representing constraints that must be met by parameters.
///
/// \tparam Type of the value stored in the constrained parameter.
template <typename T>
class ParameterConstraint {
 protected:
  ParameterConstraint() = default;

  ParameterConstraint(const ParameterConstraint &other) noexcept = default;
  ParameterConstraint(ParameterConstraint &&other) noexcept = default;
  ParameterConstraint &operator=(const ParameterConstraint &other) noexcept = default;
  ParameterConstraint &operator=(ParameterConstraint &&other) noexcept = default;

 public:
  virtual ~ParameterConstraint() = default;

  /// \brief Throw an exception if the value \p value does not satisfy the constraint.
  ///
  /// \param path
  ///   Location of the constrained parameter in the configuration tree.
  /// \param value
  ///   Value to be checked.
  virtual void checkValue(const std::string &path, const T &value) const = 0;

  /// \brief Return a JSON schema describing the constraint.
  ///
  /// The result will be merged into the existing schema used to validate the parameter.
  virtual PropertyJsonSchema jsonSchema() const = 0;
};

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_PARAMETERCONSTRAINT_H_
