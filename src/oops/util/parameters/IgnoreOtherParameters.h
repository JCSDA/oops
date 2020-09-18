/*
 * (C) Copyright 2020 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_IGNOREOTHERPARAMETERS_H_
#define OOPS_UTIL_PARAMETERS_IGNOREOTHERPARAMETERS_H_

#include "oops/util/parameters/ParameterBase.h"

namespace oops {

/// \brief Acts as a sink for unrecognized top-level configuration tree entries.
///
/// The `Parameters::validate(const eckit::Configuration &config)` function normally treats any
/// top-level configuration key whose name doesn't match the name of any registered parameter as an
/// error. Adding a member of type IgnoreOtherParameters to a Parameters subclass suppresses this
/// behavior: such unrecognized keys are simply ignored.
class IgnoreOtherParameters : public ParameterBase {
 public:
  explicit IgnoreOtherParameters(Parameters *parent = nullptr) : ParameterBase(parent) {}

  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override {}

  void serialize(eckit::LocalConfiguration &config) const override {}

  ObjectJsonSchema jsonSchema() const override;
};

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_IGNOREOTHERPARAMETERS_H_
