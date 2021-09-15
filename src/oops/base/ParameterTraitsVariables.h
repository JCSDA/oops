/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_PARAMETERTRAITSVARIABLES_H_
#define OOPS_BASE_PARAMETERTRAITSVARIABLES_H_

#include <string>
#include <vector>

#include "oops/base/Variables.h"
#include "oops/util/CompositePath.h"
#include "oops/util/parameters/ParameterTraits.h"
#include "oops/util/stringFunctions.h"

/// \file ParameterTraitsVariables.h
/// This file needs to be included before any uses of (Required/Optional)Parameter<Variables>.

namespace oops {

/// \brief Returns the list of "base" variable names (i.e. names without channel suffixes)
/// stored in \p variables.
///
/// Throws an exception if some variables have different channel suffixes than others.
std::vector<std::string> getVariableNamesWithoutChannelSuffix(const Variables &variables);

/// \brief Specialization of ParameterTraits needed for serialization and deserialization of
/// instances of Variables to/from Parameter objects.
template <>
struct ParameterTraits<Variables> {
  static boost::optional<Variables> get(util::CompositePath &path,
                                        const eckit::Configuration &config,
                                        const std::string& name) {
    if (config.has(name)) {
      return Variables(config, name);
    } else {
      return boost::none;
    }
  }

  static void set(eckit::LocalConfiguration &config,
                  const std::string &name,
                  const Variables &value) {
    // Try to remove channel indices from variable names before storing the latter in the
    // LocalConfiguration object.
    //
    // This will work if 'value' hasn't been modified after construction, but if it has, it may
    // not. For example, a sum of two sets of variables, each with a different set of channels,
    // cannot be represented by a single Configuration object. In these cases
    // getVariableNamesWithoutChannelSuffix() will throw an exception.
    const std::vector<std::string> baseVariableNames = getVariableNamesWithoutChannelSuffix(value);

    config.set(name, baseVariableNames);
    const std::vector<int> &channels = value.channels();
    if (!channels.empty()) {
      const std::string channelsAsString = util::stringfunctions::join(
            ",", channels.begin(), channels.end(), [](int n) { return std::to_string(n); });
      config.set("channels", channelsAsString);
    }
  }

  static ObjectJsonSchema jsonSchema(const std::string &name) {
    ObjectJsonSchema nameSchema = ParameterTraits<std::vector<std::string>>::jsonSchema("");
    return ObjectJsonSchema({{name, nameSchema.properties().at("")},
                             {"channels", {{"type", "[\"string\", \"integer\"]"}}}});
  }

  static std::string valueAsJson(const Variables &value);
};

}  // namespace oops

#endif  // OOPS_BASE_PARAMETERTRAITSVARIABLES_H_
