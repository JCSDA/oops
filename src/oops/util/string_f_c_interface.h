/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_STRING_F_C_INTERFACE_H_
#define OOPS_UTIL_STRING_F_C_INTERFACE_H_

#include <string>
#include <vector>

// ------------------------------------------------------------------------------
// These functions provide tools for interfacing Fortran and C++ string objects
// ------------------------------------------------------------------------------

namespace oops {

extern "C" {
  void push_string_to_vector_f(std::vector<std::string> &, const char *);
}

}  // namespace oops

#endif  // OOPS_UTIL_STRING_F_C_INTERFACE_H_
