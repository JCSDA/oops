/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_WILDCARD_H_
#define OOPS_UTIL_WILDCARD_H_

#include <string>

namespace util {

/// Returns true if \p string matches \p pattern, false otherwise.
///
/// \p pattern may contain wildcards: ? stands for any character and * stands for any number of
/// characters.
bool matchesWildcardPattern(const std::string &string, const std::string &pattern);

}  // namespace util

#endif  // OOPS_UTIL_WILDCARD_H_
