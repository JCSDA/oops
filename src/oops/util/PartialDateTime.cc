/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/dateFunctions.h"
#include "oops/util/PartialDateTime.h"


namespace util {


// -----------------------------------------------------------------------------
PartialDateTime::PartialDateTime(int year, int month, int day,
                                 int hour, int minute, int second)
    : partialdt_({year, month, day, hour, minute, second}) {
  populatedvalues_ = {year != unset_, month != unset_, day != unset_,
                      hour != unset_, minute != unset_, second != unset_};
}


// -----------------------------------------------------------------------------
PartialDateTime::PartialDateTime(std::string const &datetime_string) {
  int year, month, day, hour, minute, second;
  util::datefunctions::stringToYYYYMMDDhhmmss(datetime_string, year, month, day,
                                              hour, minute, second);
  partialdt_ = {year, month, day, hour, minute, second};
  populatedvalues_ = {year != unset_, month != unset_, day != unset_,
                      hour != unset_, minute != unset_, second != unset_};
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
bool PartialDateTime::operator==(const DateTime &other) const {
  int year, month, day, hour, minute, second;
  other.toYYYYMMDDhhmmss(year, month, day, hour, minute, second);
  std::array<int, 6> cmpdt {year, month, day, hour, minute, second};
  for (size_t i=0; i < partialdt_.size(); i++) {
    if ((populatedvalues_[i]) && (partialdt_[i] != cmpdt[i])) {
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
    if ((populatedvalues_[i]) && (partialdt_[i] != cmpdt[i])) {
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
    if ((populatedvalues_[i]) && (partialdt_[i] != cmpdt[i])) {
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
