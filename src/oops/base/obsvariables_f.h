/*
 * (C) Copyright 2024- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <stddef.h>

// ------------------------------------------------------------------------------
// These functions provide tools for interfacing Fortran and C++ string objects
// ------------------------------------------------------------------------------

namespace oops {

class ObsVariables;

extern "C" {
  oops::ObsVariables* obsvariables_empty_ctor_f();
  void obsvariables_destruct_f(oops::ObsVariables *);
  void obsvariables_push_back_f(oops::ObsVariables &, const char *);
  void obsvariables_clear_f(oops::ObsVariables &);
  size_t obsvariables_size_f(const oops::ObsVariables &);
  void obsvariables_getvariablelength_f(const oops::ObsVariables &, const size_t &, size_t &);
  void obsvariables_getvariable_f(const oops::ObsVariables &, const size_t &, size_t &,
                               const size_t &, char *);
  bool obsvariables_has_f(const oops::ObsVariables & vars, const char *);
  int obsvariables_find_f(const oops::ObsVariables & vars, const char *);
}

}  // namespace oops
