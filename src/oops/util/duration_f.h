/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_UTIL_DURATION_F_H_
#define OOPS_UTIL_DURATION_F_H_

#include "oops/util/Duration.h"

namespace util {

// -----------------------------------------------------------------------------
// These functions provide a Fortran-callable interface to Duration
// -----------------------------------------------------------------------------

extern "C" {
  util::Duration * duration_construct_f(const char str[]);
  void duration_destruct_f(util::Duration *);
  int64_t duration_int_f(const util::Duration *);
  int64_t duration_str_int_f(const char str[]);
  void duration_int_str_f(const int64_t &, char str[21]);
}

}  // namespace util
#endif  // OOPS_UTIL_DURATION_F_H_
