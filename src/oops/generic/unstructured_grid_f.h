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
  void get_nlevs_f90(int &, int &);
  void get_ncols_f90(int &, int &);
  void get_lats_f90(int &, const int &, double *);
  void get_lons_f90(int &, const int &, double *);
  void get_levs_f90(int &, const int &, double *);
  void get_cmask_f90(int &, const int &, const int &, int *);
}
}  // namespace oops

#endif  // OOPS_GENERIC_UNSTRUCTURED_GRID_F_H_
