/*
 * (C) Crown copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_RANGE_H_
#define OOPS_UTIL_RANGE_H_

#include <ostream>

namespace util {

/// \brief The range of numbers from `begin` up to but not including `end`.
template <typename T>
struct Range {
  T begin;
  T end;
};

template <typename T>
bool operator==(const Range<T> &a, const Range<T> &b) {
  return a.begin == b.begin && a.end == b.end;
}

template <typename T>
std::ostream & operator<<(std::ostream &os, const Range<T> &range) {
  os << '[' << range.begin << ", " << range.end << ')';
  return os;
}

}  // namespace util

#endif  // OOPS_UTIL_RANGE_H_
