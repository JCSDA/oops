/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_ASSOCIATIVECONTAINERS_H_
#define OOPS_UTIL_ASSOCIATIVECONTAINERS_H_

#include <map>
#include <set>
#include <utility>
#include <vector>

namespace oops {
  /// \brief Return true if \p map contains the key \p key, false otherwise.
  template<typename Key, typename Value, typename Compare, typename Alloc, typename K>
  bool contains(const std::map<Key, Value, Compare, Alloc> & map, const K & key) {
    return map.find(key) != map.end();
  }

  /// \brief Return true if \p set contains the key \p key, false otherwise.
  template<typename Key, typename Compare, typename Alloc, typename K>
  bool contains(const std::set<Key, Compare, Alloc> & set, const K & elem) {
    return set.find(elem) != set.end();
  }

  /// \brief Return the vector of keys in \p map.
  template<typename Key, typename Value, typename Compare, typename Alloc>
  std::vector<Key> keys(const std::map<Key, Value, Compare, Alloc> & map) {
    std::vector<Key> result;
    result.reserve(map.size());
    for (const std::pair<const Key, Value> &keyValue : map)
      result.push_back(keyValue.first);
    return result;
  }
}  // namespace oops

#endif  // OOPS_UTIL_ASSOCIATIVECONTAINERS_H_
