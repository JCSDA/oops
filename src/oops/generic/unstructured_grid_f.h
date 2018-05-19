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
  void get_size_f90(int &, const int &, int &);
  void get_lon_f90(int &, const int &, double *);
  void get_lat_f90(int &, const int &, double *);
  void get_area_f90(int &, const int &, double *);
  void get_vunit_f90(int &, const int &, double *);
  void get_imask_f90(int &, const int &, const int &, int *);
  void get_data_f90(int &, const int &, double *);
}
}  // namespace oops

#endif  // OOPS_GENERIC_UNSTRUCTURED_GRID_F_H_
