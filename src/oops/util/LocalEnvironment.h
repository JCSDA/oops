/*
 * (C) Copyright 2021 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_LOCALENVIRONMENT_H_
#define OOPS_UTIL_LOCALENVIRONMENT_H_

#include <map>
#include <set>
#include <string>
#include <utility>

namespace util {

/// \brief Change environment variables, restoring their original values on destruction.
class LocalEnvironment {
 public:
  /// \brief Constructor.
  LocalEnvironment() = default;

  LocalEnvironment(const LocalEnvironment &) = delete;
  LocalEnvironment(LocalEnvironment &&) = delete;
  LocalEnvironment& operator=(const LocalEnvironment &) = delete;
  LocalEnvironment& operator=(LocalEnvironment &&) = delete;

  /// \brief Restore all environment variables changed by the set() function
  /// to their original values.
  ~LocalEnvironment();

  /// Set the environment variable \p variableName to \p value.
  void set(const std::string &variableName, const std::string &value);

 private:
  std::map<std::string, std::string> variableNamesAndOriginalValues_;
  std::set<std::string> originallyUnsetVariables_;
};

}  // namespace util

#endif  // OOPS_UTIL_LOCALENVIRONMENT_H_
