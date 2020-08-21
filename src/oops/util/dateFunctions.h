/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_UTIL_DATEFUNCTIONS_H_
#define OOPS_UTIL_DATEFUNCTIONS_H_

#include <stdint.h>
#include <string>


namespace util {
  /// Namespace for non-member date-handling functions
  namespace datefunctions {

    uint64_t dateToJulian(const int year, const int month, const int day);

    void julianToDate(const uint64_t julian, int & yy, int & mm, int & dd);

    int hmsToSeconds(const int hour, const int minute, const int second);

    void secondToHms(const int seconds, int & hh, int & mm, int & ss);

    bool isLeapYear(const int year);

    bool validHhmmss(const int hour, const int minute, const int second);

    bool validYYYYMMDD(const int year, const int month, const int day);

    void stringToYYYYMMDDhhmmss(const std::string & str,
                                int & year, int & month, int & day,
                                int & hour, int & minute, int & second);

  }  // namespace datefunctions
}  // namespace util

#endif  // OOPS_UTIL_DATEFUNCTIONS_H_
