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
  void create_oobump_f90(int &, const int &, const eckit::Configuration * const *, const int &,
                         const int &, const int &, const int &);
  void delete_oobump_f90(const int &);
  void add_oobump_member_f90(const int &, const int &, const int &, const int &);
  void run_oobump_drivers_f90(const int &);
  void multiply_oobump_vbal_f90(const int &, const int &);
  void multiply_oobump_vbal_inv_f90(const int &, const int &);
  void multiply_oobump_vbal_ad_f90(const int &, const int &);
  void multiply_oobump_vbal_inv_ad_f90(const int &, const int &);
  void multiply_oobump_nicas_f90(const int &, const int &);
  void get_oobump_cv_size_f90(const int &, int &);
  void multiply_oobump_nicas_sqrt_f90(const int &, const double *, const int &);
  void multiply_oobump_nicas_sqrt_ad_f90(const int &, const int &, const double *);
  void get_oobump_param_f90(const int &, const int &, const char *, const int &);
  void set_oobump_param_f90(const int &, const int &, const char *, const int &);
}

namespace bump {
  int readEnsMember(1);
  int readPseudoEnsMember(2);
}
}  // namespace oops

#endif  // OOPS_GENERIC_OOBUMP_F_H_
