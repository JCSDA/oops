/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_UTIL_FORMATS_H_
#define OOPS_UTIL_FORMATS_H_

#include <iomanip>
#include <limits>
#include <sstream>
#include <string>

namespace util {

// -----------------------------------------------------------------------------
/// Formatting numerical types with full precision for output

template<class T>
std::string full_precision(const T & zz) {
  std::stringstream sss;
  sss << std::setprecision(std::numeric_limits<T>::digits10+1) << zz;
  return sss.str();
}

// -----------------------------------------------------------------------------

}  // namespace util

#endif  // OOPS_UTIL_FORMATS_H_
