/*
 * (C) Copyright 2018-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_MISSINGVALUES_H_
#define OOPS_UTIL_MISSINGVALUES_H_


namespace util {

// Special missing value for type T. Preferred interface, accepting T as a template parameter.
template <typename T>
const T & missingValue();

// Special missing value for type T. Backwards-compatibility interface, where the type T is
// inferred from a dummy function argument.
template <typename T>
const T & missingValue(const T & /*dummy*/) {
  return missingValue<T>();
}

}  // namespace util

#endif  // OOPS_UTIL_MISSINGVALUES_H_
