/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_UTIL_DATETIME_F_H_
#define OOPS_UTIL_DATETIME_F_H_

#include <vector>
#include "oops/util/DateTime.h"

namespace util {

// -----------------------------------------------------------------------------
// These functions provide a Fortran-callable interface to DateTime
// -----------------------------------------------------------------------------

extern "C" {
  util::DateTime * datetime_construct_f(const char str[]);
  void datetime_set_f(const char str[], util::DateTime *);
  void datetime_destruct_f(util::DateTime *);
  void datetime_string_f(const util::DateTime *, char str[21]);
  void datetime_string_io_f(const util::DateTime *, char str[17]);
  void datetime_to_yyyymmddhhmmss_f(const util::DateTime *dt,
                                    int & year, int & month, int & day,
                                    int & hour, int & minute, int & second);
  int64_t datetime_seconds_since_jan1_f(const util::DateTime *dt);
  void datetime_getints_f(const util::DateTime *, int64_t &, int &);
  void datetime_setints_f(util::DateTime *, const int64_t &, const int &);
  int64_t datetime_diff_f(const util::DateTime *, const util::DateTime *);
  void datetime_update_f(util::DateTime *, const int64_t &);
  void push_to_datetime_vector_f(std::vector<util::DateTime> *,
                                 util::DateTime *);
}

}  // namespace util
#endif  // OOPS_UTIL_DATETIME_F_H_
