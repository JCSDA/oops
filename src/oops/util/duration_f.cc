/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "oops/util/duration_f.h"

#include <cstring>
#include <string>

#include "eckit/exception/Exceptions.h"
#include "oops/util/Duration.h"

namespace util {

// -----------------------------------------------------------------------------

util::Duration * duration_construct_f(const char str[]) {
  std::string s(str);
  return new util::Duration(s);
}

// -----------------------------------------------------------------------------

void duration_destruct_f(util::Duration * dur) {
  delete dur;
}

// -----------------------------------------------------------------------------

int64_t duration_str_int_f(const char str[]) {
  std::string s(str);
  util::Duration d(s);
  return d.toSeconds();
}

// -----------------------------------------------------------------------------

int64_t duration_int_f(const util::Duration * dt) {
  return dt->toSeconds();
}

// -----------------------------------------------------------------------------

void duration_int_str_f(const int64_t & dt, char str[21]) {
  util::Duration d(dt);
  std::string s(d.toString());
  ASSERT(s.length() < 21);
//  snprintf(str, sizeof(str), s.c_str());
  std::strcpy(str, s.c_str());
}

// -----------------------------------------------------------------------------

}  // namespace util
