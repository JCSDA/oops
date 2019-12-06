/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_PARAMETERBASE_H_
#define OOPS_UTIL_PARAMETERS_PARAMETERBASE_H_

namespace eckit {
  class Configuration;
}

namespace oops {

class Parameters;

/// \brief Abstract interface of parameters that can be loaded from Configuration objects.
class ParameterBase {
 public:
  /// \brief Registers the newly created parameter in \p parent.
  explicit ParameterBase(Parameters *parent = nullptr);

  virtual ~ParameterBase() {}

  /// \brief Load the parameter's value from \p config, if present.
  virtual void deserialize(const eckit::Configuration &config) = 0;
};

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_PARAMETERBASE_H_
