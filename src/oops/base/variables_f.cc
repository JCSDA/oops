/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <string.h>
#include <string>

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/base/variables_f.h"

namespace oops {

// -----------------------------------------------------------------------------
void variables_push_back_f(oops::Variables & vars, const char * vname) {
  vars.push_back(std::string(vname));
}

// -----------------------------------------------------------------------------
size_t variables_size_f(const oops::Variables & vars) {
  return vars.size();
}

// -----------------------------------------------------------------------------
void variables_getvariable_f(const oops::Variables & vars, const size_t & jj,
                             const size_t & lcvarname, char * cvarname) {
  std::string varname = vars[jj];
  ASSERT(varname.size() < lcvarname);
  strncpy(cvarname, varname.c_str(), varname.size());
  cvarname[varname.size()] = '\0';
}

// -----------------------------------------------------------------------------

}  // namespace oops
