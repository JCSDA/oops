/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_MISSINGVALUES_H_
#define OOPS_UTIL_MISSINGVALUES_H_

#include <cstdint>
#include <string>

namespace util {
class DateTime;

const float &        missingValue(const float &);
const double &       missingValue(const double &);
const int16_t &      missingValue(const int16_t &);
const int32_t &      missingValue(const int32_t &);
const int64_t &      missingValue(const int64_t &);
const DateTime &     missingValue(const DateTime &);
const std::string &  missingValue(const std::string &);

}  // namespace util

#endif  // OOPS_UTIL_MISSINGVALUES_H_
