/*
 * (C) Copyright 2009-2020 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_UTIL_DATETIME_H_
#define OOPS_UTIL_DATETIME_H_

#include <stdint.h>
#include <iostream>
#include <string>
#include <vector>

#include "oops/util/Serializable.h"

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

class DateTime : public util::Serializable {
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

  /// sets the date given two integers: YYYYMMDD, hhmmss
  DateTime(const int YYYYMMDD, const int hhmmss);

// -- Destructor
  ~DateTime();

// -- Methods

  /// Convert the datetime to ISO 8601 format: YYYY-MM-DDThh:mm:ssZ
  std::string toString() const;

  /// Convert the datetime to format: YYYYMMDDThhmmssZ
  std::string toStringIO() const;

  // Convert the datetime to set of integers, year, month, day, hour, minute, second
  void toYYYYMMDDhhmmss(int & year, int & month, int & day,
                        int & hour, int & minute, int & second) const;

  // Convert the datetime to two integers: YYYYMMDD and hhmmss
  void toYYYYMMDDhhmmss(int & YYYYMMDD, int & hhmmss) const;

  // Convert the datetime to seconds since Jan 1 of the year
  //
  // Performance note: this method is not optimized for large numbers of calls in the same year,
  // because each invocation computes anew the DateTime for Jan 1 00:00:00 of the year. This should
  // not be a problem for computing secondsSinceJan1 O(once per obs space), but could become a
  // problem for applications that compute secondsSinceJan1 O(once per obs). In that use case, it
  // may be helpful to precompute the Jan 1 DateTime.
  int64_t secondsSinceJan1() const;

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

  bool isSet() const;

  // Serialize and deserialize
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

 private:
// -- Copy allowed
// DateTime(const DateTime&); -- default shallow copy is OK
// DateTime& operator=(const DateTime&); -- default assignment is OK

/// Set the time from an ISO 8601 format string: ${date}T${time}Z
/// where date is YYYYMMDD or YYYY-MM-DD and time is hhmmss or hh:mm:ss
  void set(const std::string &);

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

#endif  // OOPS_UTIL_DATETIME_H_
