/*
 * (C) Copyright 2019 Met Office UK
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_PARAMETERS_PARAMETERS_H_
#define OOPS_UTIL_PARAMETERS_PARAMETERS_H_

#include <set>
#include <string>
#include <vector>

#include "oops/util/parameters/ParameterBase.h"

namespace eckit {
  class Configuration;
}

namespace oops {

/// \brief A collection of parameters that can be loaded together from a \p Configuration object.
///
/// Individual parameters are represented by \p Parameter, \p OptionalParameter and
/// \p RequiredParameter objects. To register a parameter as a member of the collection
/// represented by an instance \p Parameters, pass a pointer to that instance of the parameter's
/// constructor.
class Parameters {
 public:
  /// When \p spellcheck is set to true, warning messages will be printed to
  /// \p eckit::Log::warning() if any keys in the \p Configuration object passed to
  /// \p deserialize() do not correspond to names of previously registered parameters or
  /// external parameters.
  explicit Parameters(bool spellcheck = false) : spellcheckEnabled_(spellcheck) {}

  /// \brief Add \p parameter to the list of parameters processed by subsequent calls to
  /// deserialize().
  void registerChild(ParameterBase &parameter);

  /// \brief Add \p externalParam to the list of parameter names accepted by the spellchecker.
  ///
  /// Some subclasses may encapsulate only a subset of options located on a particular level of the
  /// hierarchy represented by an \p eckit::Configuration object, with other options being read
  /// directly from that object. To prevent the names of such ("external") options from being
  /// flagged as typos by the spellchecker, they should be registered by calling
  /// \p registerExternalParameter() prior to \p deserialize().

  void registerExternalParameter(std::string externalParam);

  /// \brief Load all previously registered parameters from \p config
  void deserialize(const eckit::Configuration &config);

 private:
  std::vector<ParameterBase*> children_;
  std::set<std::string> externalParameters_;
  bool spellcheckEnabled_;
};

}  // namespace oops

#endif  // OOPS_UTIL_PARAMETERS_PARAMETERS_H_
