/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include "oops/util/PartialDateTime.h"

#include <iomanip>

#include "eckit/exception/Exceptions.h"
#include "oops/util/dateFunctions.h"


namespace util {


// -----------------------------------------------------------------------------
PartialDateTime::PartialDateTime(int year, int month, int day,
                                 int hour, int minute, int second)
    : partialdt_({year, month, day, hour, minute, second}) {

  datetime_string_ = padNumber(year, 4);
  datetime_string_ += "-" + padNumber(month, 2);
  datetime_string_ += "-" + padNumber(day, 2);
  datetime_string_ += "T" + padNumber(hour, 2);
  datetime_string_ += ":" + padNumber(minute, 2);
  datetime_string_ += ":" + padNumber(second, 2) + "Z";
}


// -----------------------------------------------------------------------------
PartialDateTime::PartialDateTime(std::string const &datetime_string) {
  int year, month, day, hour, minute, second;
  datetime_string_ = datetime_string;

  // Derive an ordinary datetime string
  std::string alt_datetime_string = datetime_string;

  // Determine if stringUnset_ is present and for which component(s).
  if (datetime_string.size() != 20)
    throw eckit::BadParameter("Partial date-time string '" + datetime_string +
                              "' is of unexpected length");

  std::array<bool, 6> populatedvalues = {true, true, true, true, true, true};
  bool isSet;
  for (size_t ind=0; ind < datetime_string.size(); ind++) {
    isSet = datetime_string[ind] != stringUnset_;
    if (!isSet) alt_datetime_string[ind] = '0';  // Replace * with 0
    if (ind <= 3) {
      populatedvalues[0] &= isSet;
    } else if (ind >= 5 && ind <= 6) {
      populatedvalues[1] &= isSet;
    } else if (ind >= 8 && ind <= 9) {
      populatedvalues[2] &= isSet;
    } else if (ind >= 11 && ind <= 12) {
      populatedvalues[3] &= isSet;
    } else if (ind >= 14 && ind <= 15) {
      populatedvalues[4] &= isSet;
    } else if (ind >= 17 && ind <= 18) {
      populatedvalues[5] &= isSet;
    }
  }
  util::datefunctions::stringToYYYYMMDDhhmmss(alt_datetime_string, year, month, day,
                                              hour, minute, second);
  // Replace with unset values
  partialdt_ = {year, month, day, hour, minute, second};
  for (size_t ind=0; ind < partialdt_.size(); ind++) {
    if (!populatedvalues[ind]) partialdt_[ind] = intUnset_;
  }
}


// -----------------------------------------------------------------------------
int PartialDateTime::year() const {return partialdt_[0];}


// -----------------------------------------------------------------------------
int PartialDateTime::month() const {return partialdt_[1];}


// -----------------------------------------------------------------------------
int PartialDateTime::day() const {return partialdt_[2];}


// -----------------------------------------------------------------------------
int PartialDateTime::hour() const {return partialdt_[3];}


// -----------------------------------------------------------------------------
int PartialDateTime::minute() const {return partialdt_[4];}


// -----------------------------------------------------------------------------
int PartialDateTime::second() const {return partialdt_[5];}


// -----------------------------------------------------------------------------
std::string PartialDateTime::padNumber(const int num, const int width) const {
  std::ostringstream ss;
  if (num == intUnset_) {
    ss << std::setw(width) << std::setfill(stringUnset_) << stringUnset_;
  } else {
    ss << std::setw(width) << std::setfill('0') << num;
  }
  return ss.str();
}


// -----------------------------------------------------------------------------
bool PartialDateTime::operator==(const PartialDateTime &other) const {
  if (this->year() != other.year()) return false;
  if (this->month() != other.month()) return false;
  if (this->day() != other.day()) return false;
  if (this->hour() != other.hour()) return false;
  if (this->minute() != other.minute()) return false;
  if (this->second() != other.second()) return false;
  return true;
}


// -----------------------------------------------------------------------------
bool PartialDateTime::operator!=(const PartialDateTime &other) const {
  return !(*this == other);
}


// -----------------------------------------------------------------------------
bool PartialDateTime::operator==(const DateTime &other) const {
  int year, month, day, hour, minute, second;
  other.toYYYYMMDDhhmmss(year, month, day, hour, minute, second);
  std::array<int, 6> cmpdt {year, month, day, hour, minute, second};
  for (size_t i=0; i < partialdt_.size(); i++) {
    if ((partialdt_[i] != intUnset_) && (partialdt_[i] != cmpdt[i])) {
      return false;
    }
  }
  return true;
}


// -----------------------------------------------------------------------------
bool PartialDateTime::operator<(const DateTime &other) const {
  int year, month, day, hour, minute, second;
  other.toYYYYMMDDhhmmss(year, month, day, hour, minute, second);
  std::array<int, 6> cmpdt {year, month, day, hour, minute, second};
  for (size_t i=0; i < partialdt_.size(); i++) {
    if ((partialdt_[i] != intUnset_) && (partialdt_[i] != cmpdt[i])) {
      return (partialdt_[i] < cmpdt[i]);
    }
  }
  return false;
}


// -----------------------------------------------------------------------------
bool PartialDateTime::operator>(const DateTime &other) const {
  int year, month, day, hour, minute, second;
  other.toYYYYMMDDhhmmss(year, month, day, hour, minute, second);
  std::array<int, 6> cmpdt {year, month, day, hour, minute, second};
  for (size_t i=0; i < partialdt_.size(); i++) {
    if ((partialdt_[i] != intUnset_) && (partialdt_[i] != cmpdt[i])) {
      return (partialdt_[i] > cmpdt[i]);
    }
  }
  return false;
}


// -----------------------------------------------------------------------------
bool PartialDateTime::operator!=(const DateTime &other) const {
  return !(*this == other);
}

// -----------------------------------------------------------------------------
bool PartialDateTime::operator<=(const DateTime &other) const {
  return (*this < other) || (*this == other);
}


// -----------------------------------------------------------------------------
bool PartialDateTime::operator>=(const DateTime &other) const {
  return (*this > other) || (*this == other);
}


// -----------------------------------------------------------------------------
bool operator<(const util::DateTime &dt, const PartialDateTime &pdt) {
  return (pdt > dt);
}


// -----------------------------------------------------------------------------
bool operator>(const DateTime &dt, const PartialDateTime &pdt) {
  return (pdt < dt);
}


// -----------------------------------------------------------------------------
bool operator==(const DateTime &dt, const PartialDateTime &pdt) {
    return (pdt == dt);
}


// -----------------------------------------------------------------------------
bool operator!=(const DateTime &dt, const PartialDateTime &pdt) {
  return (pdt != dt);
}


// -----------------------------------------------------------------------------
bool operator<=(const DateTime &dt, const PartialDateTime &pdt) {
  return (pdt >= dt);
}


// -----------------------------------------------------------------------------
bool operator>=(const DateTime &dt, const PartialDateTime &pdt) {
  return (pdt <= dt);
}


}  // namespace util
