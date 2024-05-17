/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "oops/base/ObsVariables.h"
#include "oops/util/CompositePath.h"
#include "oops/util/parameters/ParameterTraits.h"

/// \file ParameterTraitsObsVariables.h
/// This file needs to be included before any uses of (Required/Optional)Parameter<ObsVariables>.

namespace oops {

/// \brief Returns the list of "base" variable names (i.e. names without channel suffixes)
/// stored in \p variables.
///
/// Throws an exception if some variables have different channel suffixes than others.
std::vector<std::string> getVariableNamesWithoutChannelSuffix(const ObsVariables &variables);

/// \brief Specialization of ParameterTraits needed for serialization and deserialization of
/// instances of ObsVariables to/from Parameter objects.
template <>
struct ParameterTraits<ObsVariables> {
  static boost::optional<ObsVariables> get(util::CompositePath &path,
                                        const eckit::Configuration &config,
                                        const std::string& name) {
    if (config.has(name)) {
      return ObsVariables(config, name);
    } else {
      return boost::none;
    }
  }

  static void set(eckit::LocalConfiguration &config,
                  const std::string &name,
                  const ObsVariables &value) {
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

  static std::string valueAsJson(const ObsVariables &value);
};

}  // namespace oops
