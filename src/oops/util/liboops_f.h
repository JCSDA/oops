/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_UTIL_LIBOOPS_F_H_
#define OOPS_UTIL_LIBOOPS_F_H_

#include "oops/util/LibOOPS.h"

// -----------------------------------------------------------------------------
// These functions provide a Fortran-callable interface to LibOOPS
// -----------------------------------------------------------------------------

namespace oops {

extern "C" {
  void liboops_initialise_f();
  void liboops_finalise_f();
}

}  // namespace oops

#endif  // OOPS_UTIL_LIBOOPS_F_H_
