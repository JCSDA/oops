/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_COMPOSITEPATH_H_
#define OOPS_UTIL_COMPOSITEPATH_H_

#include <string>

namespace util {

/// \brief A path composed of multiple components that can be added and removed in an exception-safe
/// manner.
///
/// To append a component to the path, create a PathComponent object. The component will be removed
/// once that object goes out of scope.
///
/// This class is useful to keep track of the path to nodes of a tree processed recursively.
class CompositePath {
 public:
  /// \brief Create a path with no components.
  ///
  /// \param separator
  ///   Character used to separate path components.
  explicit CompositePath(char separator = '/');

  /// \brief Return a string representation of the path, i.e. a concatenation of its components,
  /// each preceded by the separator specified in the constructor.
  const std::string &path() const {
    return path_;
  }

 private:
  void push(const std::string &component);

  void pop(size_t componentLength);

 private:
  char separator_;
  std::string path_;

  friend class PathComponent;
};

class PathComponent {
 public:
  /// \brief Append a new component \p newComponent to the path \p path.
  PathComponent(CompositePath &path, const std::string &newComponent);

  /// \brief Remove the component from the path to which it was added in the constructor.
  ~PathComponent();

  PathComponent(const PathComponent &) = delete;
  PathComponent& operator=(const PathComponent &) = delete;

 private:
  CompositePath &path_;
  size_t length_;
};

}  // namespace util

#endif  // OOPS_UTIL_COMPOSITEPATH_H_
