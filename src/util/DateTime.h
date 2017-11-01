/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef UTIL_DATETIME_H_
#define UTIL_DATETIME_H_

#include <iostream>
#include <string>

#include <stdint.h>

// Forward declarations
namespace util {
  class Duration;
}

namespace util {

// Forward declarations
    // ---- No forward declarations ---- //

/// This class represents time, and provides methods to manipulate time
/*!
 * DateTime is represented internally as UTC quantized to the nearest second.
 */

class DateTime {
/// << and >> can be used to output and set the time as ISO 8601 strings
  friend std::istream& operator>>(std::istream& is, DateTime& uri);
  friend std::ostream& operator<<(std::ostream&, const DateTime&);

 public:

// -- Constructors
  /// Create a DateTime initialised to an invalid value.
  DateTime();

  /// Set the date given an ISO 8601 string
  explicit DateTime(const std::string &);

  /// sets the date given YYYY,MM,DD,hh,mm,ss
  DateTime(const int &, const int &, const int &,
           const int &, const int &, const int &);

// -- Destructor
  //  ~DateTime()  -- not required. This is a simple class.

// -- Methods

  /// Convert the time to ISO 8601 format: YYYY-MM-DDThh:mm:ssZ
  std::string toString() const;
  void toYYYYMMDDhhmmss(int & year, int & month, int & day,
                        int & hour, int & minute, int & second) const;

  // Operators to add/subtract a Duration to/from a DateTime
  DateTime& operator+=(const Duration &);
  DateTime& operator-=(const Duration &);
  const DateTime operator+(const Duration &) const;
  const DateTime operator-(const Duration &) const;

  // Difference in seconds between two DateTimes
  const Duration operator-(const DateTime&) const;

  // Comparison operators
  bool operator==(const DateTime&) const;
  bool operator!=(const DateTime&) const;
  bool operator<(const DateTime&) const;
  bool operator<=(const DateTime&) const;
  bool operator>(const DateTime&) const;
  bool operator>=(const DateTime&) const;
  std::size_t timestamp() const;

 private:

// -- Copy allowed
// DateTime(const DateTime&); -- default shallow copy is OK
// DateTime& operator=(const DateTime&); -- default assignment is OK

/// Set the time from an ISO 8601 format string: ${date}T${time}Z
/// where date is YYYYMMDD or YYYY-MM-DD and time is hhmmss or hh:mm:ss
  void set(const std::string &);

  void stringToYYYYMMDDhhmmss(const std::string & str,
                              int & year, int & month, int & day,
                              int & hour, int & minute, int & second) const;
  int eatChars(std::istream &, int) const;
  void failBadFormat(const std::string&) const;

  void addSeconds(const int64_t &);
  void failIfUnset() const;

// -- Members

  uint64_t date_;
  int time_;
};

// -----------------------------------------------------------------------------

inline std::ostream& operator<<(std::ostream& output, const DateTime& t) {
    output << t.toString();
    return output;
}

// -----------------------------------------------------------------------------

inline std::istream& operator>> (std::istream& input, DateTime& t ) {
    std::string time_str;
    input >> time_str;
    t.set(time_str);
    return input;
}

// -----------------------------------------------------------------------------

/// To use a date as key in a boost::map
std::size_t hash_value(const util::DateTime&);

// -----------------------------------------------------------------------------

}  // namespace util

#endif  // UTIL_DATETIME_H_
