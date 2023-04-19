/*
 * (C) Copyright 2009-2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "oops/util/dateFunctions.h"

#include <stdint.h>
#include <limits>
#include <sstream>

#include <boost/lexical_cast.hpp>

#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"


// -----------------------------------------------------------------------------

void failBadFormat(const std::string& str) {
  std::string message = "Badly formatted date: ";
  message.append(str);
  throw eckit::BadValue(message);
}


// -----------------------------------------------------------------------------

int eatChars(std::istream & is, int nchars) {
  // consume nchars characters from the stream and interpret as an integer
  if (nchars < 0) throw eckit::BadParameter("Cannot read a negative number of characters.");
  std::string str((size_t) nchars, '\0');
  is.get(&str[0], nchars+1);  // nchars+1 because istream.get reads (count-1) chars.

  int ret = 0;
  try {
    ret = std::stoi(str);
  }
  catch (...) {
    failBadFormat(str);
  }
  return ret;
}


/// Non-member functions for manipulating date and time

namespace util {
namespace datefunctions {

// -----------------------------------------------------------------------------

uint64_t dateToJulian(const int year, const int month, const int day) {
//  Compute the Julian Day number using the following formula
//  from https://doi.org/10.1145/364096.364097
//
//  julian_day = ( 1461 * ( y + 4800 + ( m - 14 ) / 12 ) ) / 4 +
//               ( 367 * ( m - 2 - 12 * ( ( m - 14 ) / 12 ) ) ) / 12 -
//                   ( 3 * ( ( y + 4900 + ( m - 14 ) / 12 ) / 100 ) ) / 4 +
//                     d - 32075
  if (year == 0 && month == 0 && day == 0) {
    return 0;
  }

  if (!util::datefunctions::validYYYYMMDD(year, month, day)) {
    oops::Log::error() << "year=" << year << " month=" << month << " day=" << day << std::endl;
    throw eckit::BadParameter("Not a valid date");
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
// Compute the year, month, and day from the Julian Day number
// using the formula from https://doi.org/10.1145/364096.364097
  uint64_t l     = 0;
  uint64_t n     = 0;
  uint64_t i     = 0;
  uint64_t j     = 0;
  uint64_t day   = 0;
  uint64_t month = 0;
  uint64_t year  = 0;

  if (julian == 0) {
    yy = 0;
    mm = 0;
    dd = 0;
  } else {
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
      oops::Log::error() << "year=" << year << std::endl;
      throw eckit::BadParameter("Year out of range");
    }
  }
}

// -----------------------------------------------------------------------------

int hmsToSeconds(const int hour, const int minute, const int second) {
  if (hour == 0 && minute == 0 && second == 0) {
    return 0;
  }

  if (!util::datefunctions::validHhmmss(hour, minute, second)) {
    oops::Log::error() << "hour=" << hour << " minute=" << minute <<
               " second=" << second << std::endl;
    throw eckit::BadParameter("Not a valid time");
  }
  return 3600 * hour + 60 * minute + second;
}

// -----------------------------------------------------------------------------

void secondToHms(const int seconds, int & hh, int & mm, int & ss) {
  if (seconds == 0) {
    hh = 0;
    mm = 0;
    ss = 0;
  } else {
    const int secondsPerDay = 86400;

    if (seconds >= 0 && seconds <= secondsPerDay) {
      int local_sec = seconds;
      hh = local_sec / 3600;
      local_sec %= 3600;
      mm = local_sec / 60;
      local_sec %= 60;
      ss = local_sec;
    } else {
      oops::Log::error() << "seconds=" << seconds << std::endl;
      throw eckit::BadParameter("seconds out of range");
    }
  }
}

// -----------------------------------------------------------------------------

bool isLeapYear(const int year) {
  return ((year%4 == 0 && year%100 != 0) || year%400 == 0);
}

// -----------------------------------------------------------------------------

bool validHhmmss(const int hour, const int minute, const int second) {
  return (hour   >= 0 && hour   <= 23 &&
          minute >= 0 && minute <= 59 &&
          second >= 0 && second <= 59) ||
         (hour == 0 && minute == 0 && second == 0);
}

// -----------------------------------------------------------------------------

bool validYYYYMMDD(const int year, const int month, const int day) {
  if (year == 0 && month == 0 && day == 0) {
    return true;
  }

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
/*! Set year, month, day, hour, minute and second from an ISO 8601 formatted string
 *  \param str ISO 8601 formatted string
 *  \param year
 *  \param month
 *  \param day
 *  \param hour
 *  \param minute
 *  \param second
 */
void stringToYYYYMMDDhhmmss(const std::string & str,
                            int & year, int & month, int & day,
                            int & hour, int & minute, int & second) {
  std::istringstream datestream(str);

  year = eatChars(datestream, 4);

  bool dashes = (datestream.peek() == '-');
  if (dashes) datestream.get();

  month = eatChars(datestream, 2);
  if (dashes) {
    char c = datestream.get();
    if (c != '-') {failBadFormat(str);}
  }

  day = eatChars(datestream, 2);

  char c = datestream.get();
  if (c != 'T') {failBadFormat(str);}

  hour = eatChars(datestream, 2);

  bool colons = (datestream.peek() == ':');
  if (colons) datestream.get();

  minute = eatChars(datestream, 2);
  if (colons) {
    char c = datestream.get();
    if (c != ':') {failBadFormat(str);}
  }

  second = eatChars(datestream, 2);

  c = datestream.get();
  if (c != 'Z') {failBadFormat(str);}

  datestream.peek();
  if (!datestream.eof()) {failBadFormat(str);}
}

// -----------------------------------------------------------------------------

}  // namespace datefunctions
}  // namespace util
