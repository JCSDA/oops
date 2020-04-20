/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <stdio.h>
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
void variables_getvariablelength_f(const oops::Variables & vars, const size_t & jj,
                                   size_t & lcvarname) {
  std::string varname = vars[jj];
  lcvarname = varname.size();
}

// -----------------------------------------------------------------------------
void variables_getvariable_f(const oops::Variables & vars, const size_t & jj,
                             size_t & lcvarname, const size_t & lfvarname,
                             char * cvarname) {
  std::string varname = vars[jj];
  lcvarname = varname.size();
  /* lfvarname is the length of the string in Fortran, which must be allocated
     by the Fortran calling routine.  In order to pass varname to Fortran it is
     first rendered as a c character array called cvarname, which according to 
     the C++ standard, must terminate with a null character.  So, the length
     of cvarname and the corresponding Fortran string must be at least one
     character longer than the length of the string in C++, namely lcvarname,
     to avoid buffer overrun.
  */
  ASSERT(lfvarname > lcvarname);
  snprintf(cvarname, lcvarname+1, "%s", varname.c_str());
}

// -----------------------------------------------------------------------------
bool variables_has_f(const oops::Variables & vars, const char * vname) {
  return vars.has(std::string(vname));
}

}  // namespace oops
