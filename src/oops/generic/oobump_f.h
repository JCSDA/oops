/*
 * (C) Copyright 2017 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_OOBUMP_F_H_
#define OOPS_GENERIC_OOBUMP_F_H_

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
extern "C" {
  void oobump_create_f90(int &, const int &, const eckit::Configuration * const *, const int &,
                         const int &, const int &, const int &, const int &, const char *);
  void oobump_delete_f90(const int &);
  void oobump_get_colocated_f90(const int &, int &);
  void oobump_get_nts_f90(const int &, int &);
  void oobump_get_cv_size_f90(const int &, int &);
  void oobump_add_member_f90(const int &, const int &, const int &, const int &);
  void oobump_remove_member_f90(const int &, const int &, const int &, const int &);
  void oobump_run_drivers_f90(const int &);
  void oobump_multiply_vbal_f90(const int &, const int &);
  void oobump_multiply_vbal_inv_f90(const int &, const int &);
  void oobump_multiply_vbal_ad_f90(const int &, const int &);
  void oobump_multiply_vbal_inv_ad_f90(const int &, const int &);
  void oobump_multiply_nicas_f90(const int &, const int &);
  void oobump_multiply_nicas_sqrt_f90(const int &, const double *, const int &);
  void oobump_multiply_nicas_sqrt_ad_f90(const int &, const int &, const double *);
  void oobump_randomize_nicas_f90(const int &, const int &);
  void oobump_get_param_f90(const int &, const int &, const char *, const int &);
  void oobump_set_param_f90(const int &, const int &, const char *, const int &);
}

namespace bump {
  int readEnsMember(1);
  int readPseudoEnsMember(2);
}
}  // namespace oops

#endif  // OOPS_GENERIC_OOBUMP_F_H_
