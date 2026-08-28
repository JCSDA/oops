/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string>
#include <vector>

#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/Variables.h"

namespace oops {

boost::optional<Variables> ParameterTraits<Variables>::get(util::CompositePath &path,
    const eckit::Configuration &config,
    const std::string& name) {
  if (config.has(name)) {
    return Variables(config, name);
  } else {
    return boost::none;
  }
}

void ParameterTraits<Variables>::set(eckit::LocalConfiguration &config,
    const std::string &name,
    const Variables &value) {
  config.set(name, value.variables());
}

ObjectJsonSchema ParameterTraits<Variables>::jsonSchema(const std::string &name) {
  ObjectJsonSchema nameSchema = ParameterTraits<std::vector<std::string>>::jsonSchema("");
  return ObjectJsonSchema({{name, nameSchema.properties().at("")}});
}

std::string ParameterTraits<Variables>::valueAsJson(const Variables &value) {
  const std::vector<std::string> varNames = value.variables();
  if (varNames.empty()) {
    return "[]";
  }
  return "["
    + util::stringfunctions::join(
      ", ", varNames.begin(), varNames.end(), [](std::string s) { return "\"" + s + "\""; })
    + "]";
}

}  // namespace oops
