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
  static int const unset_ = 0;
  std::array<int, 6> partialdt_;
  std::array<bool, 6> populatedvalues_;

 public:
  // sets the date given YYYY,MM,DD,hh,mm,ss
  PartialDateTime(int year = unset_, int month = unset_, int day = unset_,
                  int hour = unset_, int minute = unset_, int second = unset_);

  // sets the date given a string
  explicit PartialDateTime(std::string const &datetime_string);

  // Comparison operators
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
};


}  // namespace util


#endif  // OOPS_UTIL_PARTIALDATETIME_H_
