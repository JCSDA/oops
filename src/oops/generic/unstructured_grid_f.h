/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_GENERIC_UNSTRUCTURED_GRID_F_H_
#define OOPS_GENERIC_UNSTRUCTURED_GRID_F_H_

#include "eckit/config/Configuration.h"

namespace oops {
extern "C" {
  void create_ug_f90(int &);
  void delete_ug_f90(int &);
  void setup_dirac_f90(int &, const eckit::Configuration * const *);
}
}  // namespace oops

#endif  // OOPS_GENERIC_UNSTRUCTURED_GRID_F_H_
