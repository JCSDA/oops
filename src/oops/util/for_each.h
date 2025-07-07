/*
 * (C) Copyright 2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#pragma once

#include <string>
#include <tuple>
#include <type_traits>
#include <utility>

#include "atlas/array/Array.h"
#include "atlas/array/Range.h"
#include "eckit/exception/Exceptions.h"
#include "oops/util/FunctionSpaceHelpers.h"

namespace util {

enum class ExecutionPattern {
  serial,
  parallel
};

enum class IndexRange {
  include_halo,
  exclude_halo
};

namespace details {

template <typename T>
struct is_config : std::false_type {};

template <>
struct is_config<ExecutionPattern> : std::true_type {};

template <>
struct is_config<IndexRange> : std::true_type {};

template <typename T>
inline constexpr bool is_config_v = is_config<std::decay_t<T>>::value;

template <typename T>
inline constexpr bool is_not_config_v = !is_config<std::decay_t<T>>::value;

// TODO(adamsl): remove this trait when Atlas version >= 0.43.0
template <typename T, typename = void>
struct has_halo_member : std::false_type {};

template <typename T>
struct has_halo_member<T, std::void_t<decltype(&T::halo)>> : std::true_type {};

template <typename T>
constexpr bool has_halo_member_v = has_halo_member<T>::value;

template <typename T>
constexpr bool no_halo_member_v = !has_halo_member<T>::value;

template <typename... Fields>
auto make_view_tuple(std::tuple<Fields...>& fields) {
  return std::apply(
      [](auto&... args) {
        return std::make_tuple(atlas::array::make_view<double, 2>(args)...);
      },
      fields);
}

ExecutionPattern getDefaultForEachExecutionPattern();

IndexRange getDefaultForEachIndexRange();

template <typename Compare, typename... Fields>
bool allSame(Compare comp, std::tuple<Fields...>& fields) {
  return std::apply(
    [&](auto& first, auto&... rest) -> bool {
      return (comp(first, rest) && ...);
    },
    fields);
}

template <typename... Fields>
bool fieldsHaveSameShape(std::tuple<Fields...>& fields) {
  return allSame(
    [](const auto& first, const auto& other) {
      return first.shape() == other.shape();
    },
    fields);
}

template <typename... Fields>
bool fieldsHaveSameHorizontalShape(std::tuple<Fields...>& fields) {
  return allSame(
    [](const auto& first, const auto& other) {
      return first.shape(0) == other.shape(0);
    },
    fields);
}

template <typename AtlasField, std::enable_if_t<has_halo_member_v<AtlasField>, int> = 0>
bool hasContiguousOwnedPoints(const AtlasField& f) {
  // Note: This routine is used to make optimization decisions. Currently it only
  // identifies 1 of 3 possible cases, the owned points are at the beginning, and ignores
  // the cases where the owned points are in the middle or at the end of the horizontal
  // index space.
  return f.halo().appended();
}

// TODO(adamsl): remove this specialization when Atlas version >= 0.43.0
template <typename AtlasField, std::enable_if_t<no_halo_member_v<AtlasField>, int> = 0>
bool hasContiguousOwnedPoints(const AtlasField& f) {
  return true;
}

template <typename AtlasField, std::enable_if_t<has_halo_member_v<AtlasField>, int> = 0>
std::pair<atlas::idx_t, atlas::idx_t>
getContiguousHorizontalIndexSpace(const AtlasField& f, bool include_halo) {
  if (include_halo) {
    return std::make_pair(0, f.shape(0));
  } else {
    ASSERT(hasContiguousOwnedPoints(f) && f.halo().appended());
    return std::make_pair(0, f.halo().begin());
  }
}

// TODO(adamsl): remove this specialization when Atlas version >= 0.43.0
template <typename AtlasField, std::enable_if_t<!has_halo_member_v<AtlasField>, int> = 0>
std::pair<atlas::idx_t, atlas::idx_t>
getContiguousHorizontalIndexSpace(const AtlasField& f, bool include_halo) {
  auto getOwnedSizeOrFail = [](const atlas::FunctionSpace& fspace) -> atlas::idx_t {
    if (auto sc = atlas::functionspace::StructuredColumns(fspace)) {
      return sc.sizeOwned();
    } else if (auto nc = atlas::functionspace::NodeColumns(fspace);
               nc &&
               nc.mesh().nodes().has_field("tij") &&
               atlas::functionspace::CubedSphereNodeColumns(fspace)) {
      return atlas::functionspace::CubedSphereNodeColumns(fspace).sizeOwned();
    } else {
      throw eckit::NotImplemented(
        "FunctionSpace type not supported for getOwnedSizeOrFail: " + fspace.type(), Here());
      return -1;  // unreachable, but avoids compiler warning
    }
  };

  if (include_halo) {
    return std::make_pair(0, f.shape(0));
  } else {
    return std::make_pair(0, getOwnedSizeOrFail(f.functionspace()));
  }
}

}  // namespace details

template <typename Functor, typename... Fields>
void for_each_value(
  ExecutionPattern pattern,
  IndexRange range,
  Functor&& f,
  Fields&&... fields) {
  auto fieldsTuple = std::forward_as_tuple(std::forward<Fields>(fields)...);

  ASSERT(details::fieldsHaveSameShape(fieldsTuple));

  const auto& firstField = std::get<0>(fieldsTuple);
  auto viewsTuple = details::make_view_tuple(fieldsTuple);

  if (range == IndexRange::exclude_halo && !details::hasContiguousOwnedPoints(firstField)) {
    // TODO(adamsl): use a mask to exclude halo points.
    throw eckit::NotImplemented(
        "Excluding halo values is not implemented for non-contiguous halos.", Here());
  } else {
    const bool include_halo = range == IndexRange::include_halo;
    const auto[iStart, iEnd] = details::getContiguousHorizontalIndexSpace(firstField, include_halo);
    const atlas::idx_t jMax = firstField.shape(1);
    if (pattern == ExecutionPattern::parallel) {
      #pragma omp parallel for collapse(2)
      for (atlas::idx_t i = iStart; i < iEnd; ++i) {
        for (atlas::idx_t j = 0; j < jMax; ++j) {
          std::apply([&](auto&&... views) { f(views(i, j)...); }, viewsTuple);
        }
      }
    } else if (pattern == ExecutionPattern::serial) {
      for (atlas::idx_t i = iStart; i < iEnd; ++i) {
        for (atlas::idx_t j = 0; j < jMax; ++j) {
          std::apply([&](auto&&... views) { f(views(i, j)...); }, viewsTuple);
        }
      }
    } else {
      throw eckit::BadParameter("Unknown execution pattern.", Here());
    }
  }
}

template <
  typename Functor,
  typename... Fields, typename = typename std::enable_if_t< details::is_not_config_v<Functor>>
>
void for_each_value(Functor&& f, Fields&&... fields) {
  auto pattern = details::getDefaultForEachExecutionPattern();
  auto range = details::getDefaultForEachIndexRange();
  for_each_value(pattern, range, std::forward<Functor>(f), std::forward<Fields>(fields)...);
}

template <
  typename Functor,
  typename... Fields, typename = typename std::enable_if_t< details::is_not_config_v<Functor>>
>
void for_each_value(ExecutionPattern pattern, Functor&& f, Fields&&... fields) {
  auto range = details::getDefaultForEachIndexRange();
  for_each_value(pattern, range, std::forward<Functor>(f), std::forward<Fields>(fields)...);
}

template <
  typename Functor,
  typename... Fields, typename = typename std::enable_if_t< details::is_not_config_v<Functor>>
>
void for_each_value(IndexRange range, Functor&& f, Fields&&... fields) {
  auto pattern = details::getDefaultForEachExecutionPattern();
  for_each_value(pattern, range, std::forward<Functor>(f), std::forward<Fields>(fields)...);
}

template <typename Functor, typename... Fields>
void for_each_column(
  ExecutionPattern pattern,
  IndexRange range,
  Functor&& f,
  Fields&&... fields) {
  auto fieldsTuple = std::forward_as_tuple(std::forward<Fields>(fields)...);

  ASSERT(details::fieldsHaveSameHorizontalShape(fieldsTuple));

  const auto& firstField = std::get<0>(fieldsTuple);
  auto viewsTuple = details::make_view_tuple(fieldsTuple);

  if (range == IndexRange::exclude_halo && !details::hasContiguousOwnedPoints(firstField)) {
    // TODO(adamsl): use a mask to exclude halo points.
    throw eckit::NotImplemented(
        "Excluding halo values is not implemented for non-contiguous halos.", Here());
  } else {
    const bool include_halo = range == IndexRange::include_halo;
    const auto[iStart, iEnd] = details::getContiguousHorizontalIndexSpace(firstField, include_halo);
    if (pattern == ExecutionPattern::parallel) {
      #pragma omp parallel for
      for (atlas::idx_t i = iStart; i < iEnd; ++i) {
        std::apply([&](auto&&... views) {
          f(views.slice(i, atlas::array::Range::all())...);
        }, viewsTuple);
      }
    } else if (pattern == ExecutionPattern::serial) {
      for (atlas::idx_t i = iStart; i < iEnd; ++i) {
        std::apply([&](auto&&... views) {
          f(views.slice(i, atlas::array::Range::all())...);
        }, viewsTuple);
      }
    } else {
      throw eckit::BadParameter("Unknown execution pattern.", Here());
    }
  }
}

template <
  typename Functor,
  typename... Fields, typename = typename std::enable_if_t< details::is_not_config_v<Functor>>
>
void for_each_column(Functor&& f, Fields&&... fields) {
  auto pattern = details::getDefaultForEachExecutionPattern();
  auto range = details::getDefaultForEachIndexRange();
  for_each_column(pattern, range, std::forward<Functor>(f), std::forward<Fields>(fields)...);
}

template <
  typename Functor,
  typename... Fields, typename = typename std::enable_if_t< details::is_not_config_v<Functor>>
>
void for_each_column(ExecutionPattern pattern, Functor&& f, Fields&&... fields) {
  auto range = details::getDefaultForEachIndexRange();
  for_each_column(pattern, range, std::forward<Functor>(f), std::forward<Fields>(fields)...);
}

template <
  typename Functor,
  typename... Fields, typename = typename std::enable_if_t< details::is_not_config_v<Functor>>
>
void for_each_column(IndexRange range, Functor&& f, Fields&&... fields) {
  auto pattern = details::getDefaultForEachExecutionPattern();
  for_each_column(pattern, range, std::forward<Functor>(f), std::forward<Fields>(fields)...);
}

}  // namespace util
