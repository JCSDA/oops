/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/missingValues.h"

#include <limits>
#include "DateTime.h"

namespace util {

// Use missing value that is unlikely to be a useful value and not a remarquable
// number that could be used somewhere else (ie not exactly the max value)

// -----------------------------------------------------------------------------
const float & missingValue(const float &) {
  static const float fmiss = std::numeric_limits<float>::lowest() * 0.99;
  return fmiss;
}
// -----------------------------------------------------------------------------
const double & missingValue(const double &) {
  static const float fmiss = std::numeric_limits<float>::lowest() * 0.98;
  static const double dmiss = static_cast<double>(fmiss);
  return dmiss;
}
// -----------------------------------------------------------------------------
const int16_t & missingValue(const int16_t &) {
  static const int16_t imiss = std::numeric_limits<int16_t>::lowest() + 3;
  return imiss;
}
// -----------------------------------------------------------------------------
const int32_t & missingValue(const int32_t &) {
  static const int32_t imiss = std::numeric_limits<int32_t>::lowest() + 5;
  return imiss;
}
// -----------------------------------------------------------------------------
const int64_t & missingValue(const int64_t &) {
  static const int64_t imiss = std::numeric_limits<int64_t>::lowest() + 7;
  return imiss;
}
// -----------------------------------------------------------------------------
const DateTime & missingValue(const DateTime &) {
  static const DateTime dtmiss(9996, 2, 29, 23, 58, 57);
  return dtmiss;
}
// -----------------------------------------------------------------------------

}  // namespace util
