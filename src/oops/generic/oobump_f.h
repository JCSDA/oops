/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_OOBUMP_F_H_
#define OOPS_GENERIC_OOBUMP_F_H_

#include "atlas/field/FieldSet.h"
#include "atlas/functionspace/detail/FunctionSpaceImpl.h"

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
extern "C" {
  void oobump_create_f90(int &, const eckit::mpi::Comm *,
                         atlas::functionspace::FunctionSpaceImpl *, atlas::field::FieldSetImpl *,
                         const eckit::Configuration * const *,
                         const eckit::Configuration * const *,
                         const int &, const int &, const int &, const int &);
  void oobump_delete_f90(const int &);
  void oobump_get_cv_size_f90(const int &, int &);
  void oobump_add_member_f90(const int &, atlas::field::FieldSetImpl *, const int &, const int &);
  void oobump_remove_member_f90(const int &, atlas::field::FieldSetImpl *, const int &,
                                const int &);
  void oobump_run_drivers_f90(const int &);
  void oobump_multiply_vbal_f90(const int &, atlas::field::FieldSetImpl *);
  void oobump_multiply_vbal_inv_f90(const int &, atlas::field::FieldSetImpl *);
  void oobump_multiply_vbal_ad_f90(const int &, atlas::field::FieldSetImpl *);
  void oobump_multiply_vbal_inv_ad_f90(const int &, atlas::field::FieldSetImpl *);
  void oobump_multiply_nicas_f90(const int &, atlas::field::FieldSetImpl *);
  void oobump_multiply_nicas_sqrt_f90(const int &, const double *, atlas::field::FieldSetImpl *);
  void oobump_multiply_nicas_sqrt_ad_f90(const int &, atlas::field::FieldSetImpl *, const double *);
  void oobump_randomize_nicas_f90(const int &, atlas::field::FieldSetImpl *);
  void oobump_get_param_f90(const int &, const int &, const char *, atlas::field::FieldSetImpl *);
  void oobump_set_param_f90(const int &, const int &, const char *, atlas::field::FieldSetImpl *);
}

namespace bump {
  int readEnsMember(1);
  int readPseudoEnsMember(2);
}
}  // namespace oops

#endif  // OOPS_GENERIC_OOBUMP_F_H_
