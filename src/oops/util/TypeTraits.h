/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_TYPETRAITS_H_
#define OOPS_UTIL_TYPETRAITS_H_

#include <type_traits>

namespace cpp17 {

template <typename...>
using void_t = void;

}  // namespace cpp17

namespace util {

/// A type trait whose `value` member is set to true if and only if the first type in the list of
/// template parameters (`T`) is the same as any of the following types.
template <typename T, typename Head, typename... Tail>
struct any_is_same {
  static const bool value = std::is_same<T, Head>::value || any_is_same<T, Tail...>::value;
};

/// A type trait whose `value` member is set to true if and only if T is the same type as U.
template<typename T, typename U>
struct any_is_same<T, U> {
  static const bool value = std::is_same<T, U>::value;
};

}  // namespace util

#endif  // OOPS_UTIL_TYPETRAITS_H_
