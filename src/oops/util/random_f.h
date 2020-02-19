/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 *
 */

#ifndef OOPS_UTIL_RANDOM_F_H_
#define OOPS_UTIL_RANDOM_F_H_

#include <cstddef>
#include <cstdint>

namespace util {

// -----------------------------------------------------------------------------
/*! These functions provide a Fortran-callable interface to the 
 * C++ random number generators
 */
// -----------------------------------------------------------------------------

extern "C" {
  void random_uniform_float_f(const std::size_t &, const float &, const float &,
                              std::int32_t &, const std::size_t &,
                              const std::size_t &, float *);
  void random_uniform_double_f(const std::size_t &, const double &, const double &,
                               std::int32_t &, const std::size_t &,
                               const std::size_t &, double *);
  void random_uniform_int_f(const std::size_t &, const std::int32_t &,
                            const std::int32_t &, std::int32_t &,
                            const std::size_t &, const std::size_t &, std::int32_t *);
  void random_uniform_long_f(const std::size_t &, const std::int64_t &,
                             const std::int64_t &, std::int32_t &,
                             const std::size_t &, const std::size_t &, std::int64_t *);
  void random_normal_float_f(const std::size_t &, const float &, const float &,
                             std::int32_t &, const std::size_t &, float *);
  void random_normal_double_f(const std::size_t &, const double &, const double &,
                              std::int32_t &, const std::size_t &, double *);
}

}  // namespace util

#endif  // OOPS_UTIL_RANDOM_F_H_
