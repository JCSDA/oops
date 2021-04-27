/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#include "oops/util/PartialDateTime.h"

#include <algorithm>  // std::replace
#include <iomanip>
#include <regex>

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

  // Determine suitability of string format
  std::string expression = R"(^([0-9]{4}|\*{4})-([0-9]{2}|\*{2})-([0-9]{2}|\*{2}))"
                           R"(T([0-9]{2}|\*{2}):([0-9]{2}|\*{2}):([0-9]{2}|\*{2})Z$)";
  static const std::regex regex(expression);
  std::smatch matches;
  if (!std::regex_match(datetime_string, matches, regex))
      throw eckit::BadParameter("Partial date-time string '" + datetime_string +
                                "' is not of the expected format '" + expression + "'");

  // Determine which components are missing.
  std::array<size_t, 6> component_index = {0, 5, 8, 11, 14, 17};
  std::array<bool, 6> populatedvalues = {true, true, true, true, true, true};
  for (size_t ind=0; ind < component_index.size(); ind++)
      populatedvalues[ind] = datetime_string[component_index[ind]] != stringUnset_;

  // Replace these with '0' so that we can determine the year, month, etc.
  std::replace(alt_datetime_string.begin(), alt_datetime_string.end(), char{stringUnset_}, '0');

  util::datefunctions::stringToYYYYMMDDhhmmss(alt_datetime_string, year, month, day,
                                              hour, minute, second);
  // Replace with unset values
  partialdt_ = {year, month, day, hour, minute, second};
  for (size_t ind=0; ind < partialdt_.size(); ind++)
    if (!populatedvalues[ind]) partialdt_[ind] = intUnset_;
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
