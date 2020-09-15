/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_PARAMETERBASE_H_
#define OOPS_UTIL_PARAMETERS_PARAMETERBASE_H_

#include <set>
#include <string>

namespace eckit {
  class Channel;
  class Configuration;
  class LocalConfiguration;
}

namespace util {
  class CompositePath;
}

namespace oops {

class ObjectJsonSchema;
class Parameters;

/// \brief Abstract interface of parameters that can be loaded from and saved to Configuration
/// objects.
class ParameterBase {
 protected:
  /// \brief Registers the newly created parameter in \p parent.
  explicit ParameterBase(Parameters *parent = nullptr);

  ParameterBase(const ParameterBase &other) noexcept = default;
  ParameterBase(ParameterBase &&other) noexcept = default;
  ParameterBase &operator=(const ParameterBase &other) noexcept = default;
  ParameterBase &operator=(ParameterBase &&other) noexcept = default;

 public:
  virtual ~ParameterBase() = default;

  /// \brief Load the parameter's value from \p config, if present.
  ///
  /// \param path
  ///   Location of the configuration \p config in the full configuration loaded from a YAML file.
  ///   This object may be modified in-place during deserialization, but must be restored to its
  ///   original state before this function returns.
  /// \param config
  ///   Configuration from which parameter values are to be loaded.
  virtual void deserialize(util::CompositePath &path,
                           const eckit::Configuration &config) = 0;

  /// \brief Save the parameter's value to \p config.
  virtual void serialize(eckit::LocalConfiguration &config) const = 0;

  /// \brief Return an object encapsulating the JSON schema specifying the expected structure of
  /// the JSON (or YAML) node from which this parameter is loaded.
  virtual ObjectJsonSchema jsonSchema() const = 0;
};

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_PARAMETERBASE_H_
