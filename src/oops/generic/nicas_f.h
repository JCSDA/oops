/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_GENERIC_NICAS_F_H_
#define OOPS_GENERIC_NICAS_F_H_

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
extern "C" {
  void create_nicas_f90(int &, const eckit::Configuration * const *);
  void delete_nicas_f90(const int &);
  void nicas_multiply_f90(const int &, const int &);
}
}  // namespace oops

#endif  // OOPS_GENERIC_NICAS_F_H_
