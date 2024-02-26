/*
 * (C) Copyright 2018-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/missingValues.h"

#include <cstdint>
#include <limits>
#include <string>

#include "DateTime.h"

namespace util {

// Use missing value that is unlikely to be a useful value and not a remarquable
// number that could be used somewhere else (ie not exactly the max value)

// -----------------------------------------------------------------------------
template <>
const float & missingValue<float>() {
  static const float fmiss = std::numeric_limits<float>::lowest() * 0.99;
  return fmiss;
}
// -----------------------------------------------------------------------------
template <>
const double & missingValue<double>() {
  static const float fmiss = std::numeric_limits<float>::lowest() * 0.98;
  static const double dmiss = static_cast<double>(fmiss);
  return dmiss;
}
// -----------------------------------------------------------------------------
template <>
const int16_t & missingValue<int16_t>() {
  static const int16_t imiss = std::numeric_limits<int16_t>::lowest() + 3;
  return imiss;
}
// -----------------------------------------------------------------------------
template <>
const int32_t & missingValue<int32_t>() {
  static const int32_t imiss = std::numeric_limits<int32_t>::lowest() + 5;
  return imiss;
}
// -----------------------------------------------------------------------------
template <>
const int64_t & missingValue<int64_t>() {
  static const int64_t imiss = std::numeric_limits<int64_t>::lowest() + 7;
  return imiss;
}
// -----------------------------------------------------------------------------
template <>
const bool & missingValue<bool>() {
  // No sensible missingValue option for bool
  static const bool bmiss = false;
  return bmiss;
}
// -----------------------------------------------------------------------------
template <>
const char & missingValue<char>() {
  // No sensible missingValue option for char
  static const char chmiss = '\0';
  return chmiss;
}
// -----------------------------------------------------------------------------
template <>
const std::string & missingValue<std::string>() {
  static const std::string strmiss("MISSING*");
  return strmiss;
}
// -----------------------------------------------------------------------------
template <>
const DateTime & missingValue<DateTime>() {
  static const DateTime dtmiss(9996, 2, 29, 23, 58, 57);
  return dtmiss;
}
// -----------------------------------------------------------------------------

}  // namespace util
