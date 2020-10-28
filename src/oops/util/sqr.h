/*
 * (C) Copyright 2020 UK Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_SQR_H_
#define OOPS_UTIL_SQR_H_

namespace util {

/// \brief Return the square of \p x.
template <typename T>
T sqr(T x) {
  return x * x;
}

}  // namespace util

#endif  // OOPS_UTIL_SQR_H_

