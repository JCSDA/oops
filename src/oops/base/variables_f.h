/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_VARIABLES_F_H_
#define OOPS_BASE_VARIABLES_F_H_

#include <string>

#include "oops/base/Variables.h"

// ------------------------------------------------------------------------------
// These functions provide tools for interfacing Fortran and C++ string objects
// ------------------------------------------------------------------------------

namespace oops {

extern "C" {
  void variables_push_back_f(oops::Variables &, const char *);
  size_t variables_size_f(const oops::Variables &);
  void variables_getvariablelength_f(const oops::Variables &, const size_t &, size_t &);
  void variables_getvariable_f(const oops::Variables &, const size_t &, size_t &,
                               const size_t &, char *);
}

}  // namespace oops

#endif  // OOPS_BASE_VARIABLES_F_H_
