/*
 * (C) Copyright 2020 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_CONFIGURATIONPARAMETER_H_
#define OOPS_UTIL_PARAMETERS_CONFIGURATIONPARAMETER_H_

#include <set>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/value/Value.h"
#include "oops/util/parameters/ParameterBase.h"

namespace oops {

/// \brief A parameter storing a copy of the eckit::Configuration object from which its parent
/// Parameters object was deserialized.
///
/// This can be useful when some options located at a given level of the configuration tree have
/// already been mapped to dedicated Parameter (or RequiredParameter/OptionalParameter) member
/// variables of a Parameters subclass, but others haven't. Adding a member variable of type
/// ConfigurationParameter to that subclass will then make it possible for users of that subclass
/// to access the values of options not deserialized into other member variables.
class ConfigurationParameter : public ParameterBase {
 public:
  explicit ConfigurationParameter(Parameters *parent) : ParameterBase(parent) {}

  /// \brief Stores a copy of \p config in this object.
  void deserialize(util::CompositePath &path, const eckit::Configuration &config) override;

  void serialize(eckit::LocalConfiguration &config) const override;

  ObjectJsonSchema jsonSchema() const override;

  /// \brief The value stored in this parameter.
  const eckit::LocalConfiguration &value() const { return value_; }

  /// \brief The value stored in this parameter.
  operator const eckit::LocalConfiguration &() const { return value_; }

 private:
  std::string name_;
  eckit::LocalConfiguration value_;
};

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_CONFIGURATIONPARAMETER_H_
