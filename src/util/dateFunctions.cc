/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "util/dateFunctions.h"
#include "util/abor1_cpp.h"
#include "util/Logger.h"

#include <limits>
#include <stdint.h>

/// Non-member functions for manipulating date and time

using oops::Log;

namespace util {
namespace datefunctions {

// -----------------------------------------------------------------------------

uint64_t dateToJulian(const int year, const int month, const int day) {
//  Compute the Julian Day number applying the following formula
//
//  julian_day = ( 1461 * ( y + 4800 + ( m - 14 ) / 12 ) ) / 4 +
//               ( 367 * ( m - 2 - 12 * ( ( m - 14 ) / 12 ) ) ) / 12 -
//                   ( 3 * ( ( y + 4900 + ( m - 14 ) / 12 ) / 100 ) ) / 4 +
//                     d - 32075

  if (!util::datefunctions::validYYYYMMDD(year, month, day)) {
    Log::error() << "year=" << year << " month=" << month << " day=" << day << std::endl;
    ABORT("Not a valid date");
  }

  int m1 = (month - 14)/12;
  int a  = (1461 * (year + 4800 + m1))/4;
  int b  = (367 * (month - 2 - (12 * m1)))/12;
  int m2 = (year + 4900 + m1)/100;
  int c  = (3 * (m2))/4;
  return static_cast<uint64_t>(a + b - c + day - 32075);
}

// -----------------------------------------------------------------------------

void julianToDate(const uint64_t julian, int & yy, int & mm, int & dd) {
  uint64_t l     = 0;
  uint64_t n     = 0;
  uint64_t i     = 0;
  uint64_t j     = 0;
  uint64_t day   = 0;
  uint64_t month = 0;
  uint64_t year  = 0;

// Modified Julian date

  l = julian + static_cast<uint64_t>(68569);
  n = (4 * l) / 146097;
  l = l - (146097 * n + 3) / 4;
  i = (4000 * (l + 1)) / 1461001;
  l = l - (1461 * i) / 4 + 31;
  j = (80 * l) / 2447;
  day = l - (2447 * j) / 80;
  l = j / 11;
  month = j + 2 - (12 * l);
  year = 100 * (n - 49) + i + l;
  if (year <= static_cast<uint64_t>(std::numeric_limits<int>::max())) {
    dd = static_cast<int>(day);
    mm = static_cast<int>(month);
    yy = static_cast<int>(year);
  } else {
    Log::error() << "year=" << year << std::endl;
    ABORT("year out of range");
  }
}

// -----------------------------------------------------------------------------

int hmsToSeconds(const int hour, const int minute, const int second) {
  if (!util::datefunctions::validHhmmss(hour, minute, second)) {
    Log::error() << "hour=" << hour << " minute=" << minute <<
               " second=" << second << std::endl;
    ABORT("Not a valid time");
  }
  return 3600 * hour + 60 * minute + second;
}

// -----------------------------------------------------------------------------

void secondToHms(const int seconds, int & hh, int & mm, int & ss) {
  const int secondsPerDay = 86400;

  if (seconds >= 0 && seconds <= secondsPerDay) {
    int local_sec = seconds;
    hh = local_sec / 3600;
    local_sec %= 3600;
    mm = local_sec / 60;
    local_sec %= 60;
    ss = local_sec;
  } else {
    Log::error() << "seconds=" << seconds << std::endl;
    ABORT("seconds out of range");
  }
}

// -----------------------------------------------------------------------------

bool isLeapYear(const int year) {
  return ((year%4 == 0 && year%100 != 0) || year%400) == 0;
}

// -----------------------------------------------------------------------------

bool validHhmmss(const int hour, const int minute, const int second) {
  return hour   >= 0 && hour   <= 23 &&
         minute >= 0 && minute <= 59 &&
         second >= 0 && second <= 59;
}

// -----------------------------------------------------------------------------

bool validYYYYMMDD(const int year, const int month, const int day) {
  bool good;
  good = year   >= 0 && year <= 9999 &&
         month  >= 1 && month  <= 12 &&
         day    >= 1 && day    <= 31;

  if (good) {
    if (month == 9 || month == 4 || month == 6 || month == 11) {
        good = (day <= 30);
    } else if (month != 2) {
        good = (day <= 31);
    } else if (util::datefunctions::isLeapYear(year)) {
        good = (day <= 29);
    } else {
        good = (day <= 28);
    }
  }

  return good;
}

// -----------------------------------------------------------------------------

}  // namespace datefunctions
}  // namespace util
