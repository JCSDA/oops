/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_VARIABLES_F_H_
#define OOPS_BASE_VARIABLES_F_H_

#include <stddef.h>

// ------------------------------------------------------------------------------
// These functions provide tools for interfacing Fortran and C++ string objects
// ------------------------------------------------------------------------------

namespace oops {

class Variables;

extern "C" {
  oops::Variables* variables_empty_ctor_f();
  void variables_destruct_f(oops::Variables *);
  void variables_push_back_f(oops::Variables &, const char *);
  void variables_clear_f(oops::Variables &);
  size_t variables_size_f(const oops::Variables &);
  void variables_getvariablelength_f(const oops::Variables &, const size_t &, size_t &);
  void variables_getvariable_f(const oops::Variables &, const size_t &, size_t &,
                               const size_t &, char *);
  bool variables_has_f(const oops::Variables & vars, const char *);
  int variables_find_f(const oops::Variables & vars, const char *);
}

}  // namespace oops

#endif  // OOPS_BASE_VARIABLES_F_H_
