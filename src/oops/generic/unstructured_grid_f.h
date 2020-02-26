/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_UNSTRUCTURED_GRID_F_H_
#define OOPS_GENERIC_UNSTRUCTURED_GRID_F_H_

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/detail/FunctionSpaceImpl.h"

#include "eckit/config/Configuration.h"

#include "oops/util/DateTime.h"

namespace oops {
extern "C" {
  void create_ug_f90(int &, const int &, const int &);
  void delete_ug_f90(int &);
  void ug_get_ngrid_f90(const int &, int &);
  void ug_get_dims_f90(const int &, const int &, int &, int &, int &);
  void ug_create_atlas_grid_conf_f90(const int &, const eckit::Configuration * const *);
  void ug_set_atlas_functionspace_pointer_f90(const int &,
                                              atlas::functionspace::FunctionSpaceImpl *);
  void ug_fill_atlas_fieldset_f90(const int &, atlas::field::FieldSetImpl *);
  void ug_set_atlas_fieldset_pointer_f90(const int &,
                                         atlas::field::FieldSetImpl *);
  void ug_set_atlas_f90(const int &, atlas::field::FieldSetImpl *);
  void ug_to_atlas_f90(const int &, atlas::field::FieldSetImpl *);
  void ug_from_atlas_f90(const int &, atlas::field::FieldSetImpl *);
}
}  // namespace oops

#endif  // OOPS_GENERIC_UNSTRUCTURED_GRID_F_H_
