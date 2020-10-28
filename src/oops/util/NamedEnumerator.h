/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_NAMEDENUMERATOR_H_
#define OOPS_UTIL_NAMEDENUMERATOR_H_

namespace util {

/// \brief Value of an enumerator defined in the enumeration \c EnumType accompanied by a
/// human-readable textual representation of this value.
template <typename EnumType>
struct NamedEnumerator {
  constexpr NamedEnumerator(EnumType value_, const char* name_) :
    value(value_), name(name_)
  {}

  EnumType value;
  const char *name;
};

}  // namespace util

#endif  // OOPS_UTIL_NAMEDENUMERATOR_H_
