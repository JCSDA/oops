/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_SCALARORMAP_H_
#define OOPS_UTIL_SCALARORMAP_H_

#include <map>
#include <utility>

namespace util {

/// \brief A container that can store either a scalar of type Value or an instance of
/// std::map<Key, Value>.
///
/// If a map is stored, calls to at() and operator[]() are forwarded to the stored map.
/// If a scalar S is stored, these calls behave as if the container was storing a map with entries
/// corresponding to all possible keys set to S.
template <typename Key, typename Value>
class ScalarOrMap {
 public:
  typedef ScalarOrMap<Key, Value> Self;

  typedef std::map<Key, Value> Map;

  /// Default creates a container storing an empty map.
  ScalarOrMap() = default;

  /// Creates a container storing the map \p map.
  explicit ScalarOrMap(Map map)
    : isScalar_(false), map_(std::move(map))
  {}

  /// Creates a container storing the scalar \p scalar.
  explicit ScalarOrMap(Value scalar)
    : isScalar_(true), scalar_(std::move(scalar))
  {}

  /// Assigns a map.
  Self &operator=(Map &&map) {
    map_ = std::move(map);
    isScalar_ = false;
    return *this;
  }

  /// \overload
  Self &operator=(const Map &map) {
    map_ = map;
    isScalar_ = false;
    return *this;
  }

  /// Assigns a scalar.
  Self &operator=(Value &&scalar) {
    scalar_ = std::move(scalar);
    isScalar_ = true;
    map_.clear();
    return *this;
  }

  /// \overload
  Self &operator=(const Value &scalar) {
    scalar_ = scalar;
    isScalar_ = true;
    map_.clear();
    return *this;
  }

  /// Returns true if the container is currently storing a scalar and false if it is storing a map.
  bool isScalar() const {
    return isScalar_;
  }

  /// Returns true if the container is storing a scalar, or if it is storing a map containing the
  /// key \p key.
  bool contains(const Key &key) const {
    if (isScalar_)
      return true;
    else
      return map_.find(key) != map_.end();
  }

  /// If the container is storing a scalar, returns this scalar. If it is storing a map,
  /// returns the value corresponding to the key \p key. If the map contains no such key, an
  /// exception is thrown.
  const Value &at(const Key &key) const {
    if (isScalar_)
      return scalar_;
    else
      return map_.at(key);
  }

  /// Returns the iterator to the first element of the stored map. If the container is currently
  /// storing a scalar, an exception is thrown.
  typename Map::iterator begin() {
    if (isScalar_)
      throw std::runtime_error("The container is storing a scalar, not a map");
    return map_.begin();
  }

  /// \overload
  typename Map::const_iterator begin() const {
    if (isScalar_)
      throw std::runtime_error("The container is storing a scalar, not a map");
    return map_.begin();
  }

  /// Returns the past-the-end iterator of the stored map. If the container is currently
  /// storing a scalar, an exception is thrown.
  typename Map::iterator end() {
    if (isScalar_)
      throw std::runtime_error("The container is storing a scalar, not a map");
    return map_.end();
  }

  typename Map::const_iterator end() const {
    if (isScalar_)
      throw std::runtime_error("The container is storing a scalar, not a map");
    return map_.end();
  }

 private:
  bool isScalar_ = false;
  Value scalar_ = {};
  std::map<Key, Value> map_;
};
}  // namespace util

#endif  // OOPS_UTIL_SCALARORMAP_H_
