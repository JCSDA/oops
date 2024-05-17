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

/// \file ParameterTraitsVariables.h
/// This file needs to be included before any uses of (Required/Optional)Parameter<Variables>.

namespace oops {

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
    config.set(name, value.variables());
  }

  static ObjectJsonSchema jsonSchema(const std::string &name) {
    ObjectJsonSchema nameSchema = ParameterTraits<std::vector<std::string>>::jsonSchema("");
    return ObjectJsonSchema({{name, nameSchema.properties().at("")}});
  }

  static std::string valueAsJson(const Variables &value);
};

}  // namespace oops

#endif  // OOPS_BASE_PARAMETERTRAITSVARIABLES_H_
