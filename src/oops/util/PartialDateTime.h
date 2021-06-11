/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_PARTIALDATETIME_H_
#define OOPS_UTIL_PARTIALDATETIME_H_

#include <array>
#include <string>

#include "oops/util/DateTime.h"

namespace util {


/*!
 * PartialDateTime is a class for performing partial date/time comparisons with
 * DateTime objects.
 */


class PartialDateTime {
 private:
  std::string datetime_string_;
  static char const stringUnset_ = '*';
  static int const intUnset_ = -1;
  std::array<int, 6> partialdt_;

 public:
  // \brief Sets the date given YYYY,MM,DD,hh,mm,ss
  PartialDateTime(int year = intUnset_, int month = intUnset_, int day = intUnset_,
                  int hour = intUnset_, int minute = intUnset_, int second = intUnset_);

  // \brief Sets the date given a string
  // \brief String must represent an extended ISO 8601 format, where an asterisk is interpreted
  // as those components which are to be ignored by the comparison operators.
  explicit PartialDateTime(std::string const &datetime_string);

  // \brief Generate extended ISO 8601 string representation for this PartialDateTime
  std::string toString() const {return datetime_string_;}

  // Comparison operators
  bool operator==(const PartialDateTime &) const;
  bool operator!=(const PartialDateTime &) const;
  bool operator==(const util::DateTime &) const;
  bool operator!=(const util::DateTime &) const;
  bool operator<(const util::DateTime &) const;
  bool operator<=(const util::DateTime &) const;
  bool operator>(const util::DateTime &) const;
  bool operator>=(const util::DateTime &) const;
  friend bool operator==(const util::DateTime &dt, const PartialDateTime &pdt);
  friend bool operator!=(const util::DateTime &dt, const PartialDateTime &pdt);
  friend bool operator<(const util::DateTime &dt, const PartialDateTime &pdt);
  friend bool operator<=(const util::DateTime &dt, const PartialDateTime &pdt);
  friend bool operator>(const util::DateTime &dt, const PartialDateTime &pdt);
  friend bool operator>=(const util::DateTime &dt, const PartialDateTime &pdt);

  int year() const;
  int month() const;
  int day() const;
  int hour() const;
  int minute() const;
  int second() const;

 private:
  // \brief Used for creating a PartialDateTime string
  std::string padNumber(const int num, const int width) const;
};


}  // namespace util


#endif  // OOPS_UTIL_PARTIALDATETIME_H_
