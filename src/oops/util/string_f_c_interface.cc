/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include <string>
#include <vector>

#include "oops/util/string_f_c_interface.h"

namespace oops {

// -----------------------------------------------------------------------------
  void push_string_to_vector_f(std::vector<std::string> & vec, const char * vname) {
    vec.push_back(std::string(vname));
  }

// -----------------------------------------------------------------------------

}  // namespace oops
