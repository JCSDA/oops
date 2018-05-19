/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "util/datetime_f.h"

#include <cstring>
#include <string>

#include "eckit/exception/Exceptions.h"
#include "util/dateFunctions.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/Logger.h"

using oops::Log;

namespace util {

// -----------------------------------------------------------------------------

util::DateTime * datetime_construct_f(const char str[]) {
  std::string s(str);
  return new util::DateTime(s);
}

// -----------------------------------------------------------------------------

void datetime_set_f(const char str[], util::DateTime * dd) {
  std::string s(str);
  const util::DateTime dtmp(s);
  *dd = dtmp;
}

// -----------------------------------------------------------------------------

void datetime_destruct_f(util::DateTime * dd) {
  delete dd;
}

// -----------------------------------------------------------------------------

void datetime_string_f(const util::DateTime * dd, char cstr[21]) {
  const std::string ss(dd->toString());
  ASSERT(ss.length() < 21);
//  snprintf(cstr, sizeof(cstr), ss.c_str());
  std::strcpy(cstr, ss.c_str());
}

// -----------------------------------------------------------------------------

void datetime_getints_f(const util::DateTime * dt, int64_t & date, int & time) {
  int yy, mm, dd, hh, nn, ss;
  Log::debug() << "datetime_getints dt = " << *dt << std::endl;
  dt->toYYYYMMDDhhmmss(yy, mm, dd, hh, nn, ss);
  date = 10000 * yy + 100 * mm + dd;
  time = datefunctions::hmsToSeconds(hh, nn, ss);
  Log::debug() << "datetime_getints date, time = " << date << ", " << time << std::endl;
}

// -----------------------------------------------------------------------------

void datetime_setints_f(util::DateTime * dt, const int64_t & date, const int & time) {
  int yy, mm, dd, hh, nn, ss, ii;
  Log::debug() << "datetime_setints date, time = " << date << ", " << time << std::endl;
  yy = date / 10000;
  ii = date % 10000;
  mm = ii / 100;
  dd = ii % 100;
  datefunctions::secondToHms(time, hh, nn, ss);
  *dt = util::DateTime(yy, mm, dd, hh, nn, ss);
  Log::debug() << "datetime_setints dt = " << *dt << std::endl;
}

// -----------------------------------------------------------------------------

int64_t datetime_diff_f(const util::DateTime * dt1,
                        const util::DateTime * dt2) {
  const util::Duration dd = *dt1 - *dt2;
  return dd.toSeconds();
}

// -----------------------------------------------------------------------------

void datetime_update_f(util::DateTime * tt, const int64_t & dt) {
  const util::Duration dd(dt);
  *tt += dd;
}

// -----------------------------------------------------------------------------

}  // namespace util
