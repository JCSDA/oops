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
                                        const std::string& name);

  static void set(eckit::LocalConfiguration &config,
                  const std::string &name,
                  const ObsVariables &value);

  static ObjectJsonSchema jsonSchema(const std::string &name);

  static std::string valueAsJson(const ObsVariables &value);
};

}  // namespace oops
