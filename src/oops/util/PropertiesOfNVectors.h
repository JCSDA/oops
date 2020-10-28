/*
 * (C) Crown copyright 2020, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_PROPERTIESOFNVECTORS_H_
#define OOPS_UTIL_PROPERTIESOFNVECTORS_H_

#include <string>
#include <vector>

namespace oops {

  /// Return (as a string) a comma-delimited list of vector sizes.
  /// This is the base case for one vector.
  template <typename T>
    std::string listOfVectorSizes(const std::vector <T> &vec)
    {
      return std::to_string(vec.size());
    }

  /// Return (as a string) a comma-delimited list of vector sizes.
  /// This is the recursive case that accepts an arbitrary number of vectors
  /// using a variadic template.
  template <typename T, typename... Args>
    std::string listOfVectorSizes(const std::vector <T> &vec1,
                                  const Args&... vecs)
  {
    return std::to_string(vec1.size()) + ", " + listOfVectorSizes(vecs...);
  }
}  // namespace oops
#endif  // OOPS_UTIL_PROPERTIESOFNVECTORS_H_
