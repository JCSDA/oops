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
