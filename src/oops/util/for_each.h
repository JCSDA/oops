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
  const bool all_owned = (f.halo().size() == 0);
  const bool halo_appended = f.halo().appended();  // returns false if field has 0 halo points
  return all_owned || halo_appended;
}

// TODO(adamsl): remove this specialization when Atlas version >= 0.43.0
template <typename AtlasField, std::enable_if_t<no_halo_member_v<AtlasField>, int> = 0>
bool hasContiguousOwnedPoints(const AtlasField& f) {
  // Note: This routine is used to make optimization decisions. Currently it assumes that
  // all StructureColumn and CubedSphereNodeColumns function spaces have contiguous owned
  // points.
  const auto& fspace = f.functionspace();
  if (auto sc = atlas::functionspace::StructuredColumns(fspace)) {
    return true;
  } else if (auto nc = atlas::functionspace::NodeColumns(fspace);
              nc &&
              nc.mesh().nodes().has_field("tij") &&
              atlas::functionspace::CubedSphereNodeColumns(fspace)) {
    return true;
  } else {
    return false;
  }
}

template <typename AtlasField, std::enable_if_t<has_halo_member_v<AtlasField>, int> = 0>
std::pair<atlas::idx_t, atlas::idx_t>
getContiguousHorizontalIndexSpace(const AtlasField& f, bool include_halo) {
  if (include_halo) {
    return std::make_pair(0, f.shape(0));
  } else {
    ASSERT(hasContiguousOwnedPoints(f));
    const bool all_owned = (f.halo().size() == 0);
    if (all_owned) {
      return std::make_pair(0, f.shape(0));
    } else {
      ASSERT(f.halo().appended());  // returns false if field has 0 halo points
      return std::make_pair(0, f.halo().begin());
    }
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
        "FunctionSpace type not assumed to have contiguous owned points: " + fspace.type(), Here());
    }
  };

  if (include_halo) {
    return std::make_pair(0, f.shape(0));
  } else {
    ASSERT(hasContiguousOwnedPoints(f));

    return std::make_pair(0, getOwnedSizeOrFail(f.functionspace()));
  }
}

class ComputeSpectralCoefficientIndex {
 public:
  ComputeSpectralCoefficientIndex(
    atlas::array::LocalView<const int, 1> zonal_wavenumbers,
    int truncation);

  std::pair<atlas::idx_t, atlas::idx_t> operator()(int n, int jm, int m);
 private:
  atlas::array::ArrayT<int> prefix_sum;
  atlas::array::ArrayView<int, 1> prefix_sum_view;
};

}  // namespace details

struct IndexSpace1D {
  atlas::idx_t start;
  atlas::idx_t end;
};

struct IndexSpace2D {
  IndexSpace1D i;
  IndexSpace1D j;

  IndexSpace2D(atlas::idx_t i_start, atlas::idx_t i_end,
               atlas::idx_t j_start, atlas::idx_t j_end)
    : i{i_start, i_end}, j{j_start, j_end} {}

  IndexSpace2D(IndexSpace1D i_range, IndexSpace1D j_range)
    : i(i_range), j(j_range) {}
};

struct IndexSpace3D {
  IndexSpace1D i;
  IndexSpace1D j;
  IndexSpace1D k;

  IndexSpace3D(atlas::idx_t i_start, atlas::idx_t i_end,
               atlas::idx_t j_start, atlas::idx_t j_end,
               atlas::idx_t k_start, atlas::idx_t k_end)
    : i{i_start, i_end}, j{j_start, j_end}, k{k_start, k_end} {}

  IndexSpace3D(IndexSpace1D i_range, IndexSpace1D j_range, IndexSpace1D k_range)
    : i(i_range), j(j_range), k(k_range) {}
};

class IndexSpaceSpectral1D {
 public:
  explicit IndexSpaceSpectral1D(const atlas::functionspace::Spectral& sp)
    : zonal_wavenumbers_(sp.zonal_wavenumbers()),
      truncation_(sp.truncation()) {}

  atlas::array::LocalView<const int, 1> zonal_wavenumbers() const { return zonal_wavenumbers_; }
  int zonal_wavenumbers(int jm) const { return zonal_wavenumbers_(jm); }
  int truncation() const { return truncation_; }

 private:
  atlas::array::LocalView<const int, 1> zonal_wavenumbers_;
  int truncation_;
};

IndexSpace1D make_index_space_1d(IndexRange range, const atlas::Field& field);

IndexSpace2D make_index_space_2d(IndexRange range, const atlas::Field& field);

IndexSpace3D make_index_space_3d(IndexRange range, const atlas::Field& field);

IndexSpaceSpectral1D make_index_space_spectral_1d(const atlas::Field& field);

template <typename Functor>
void for_each_index(
  ExecutionPattern pattern,
  IndexSpace1D range,
  Functor&& f) {
  if (pattern == ExecutionPattern::parallel) {
      // Disable OpenMP parallelization for old Intel compilers to avoid compiler errors
#ifndef __INTEL_COMPILER
    #pragma omp parallel for
#endif
    for (atlas::idx_t i = range.start; i < range.end; ++i) {
      f(i);
    }
  } else if (pattern == ExecutionPattern::serial) {
    for (atlas::idx_t i = range.start; i < range.end; ++i) {
      f(i);
    }
  } else {
    throw eckit::BadParameter("Unknown execution pattern.", Here());
  }
}

template <typename Functor>
void for_each_index(
  IndexSpace1D range,
  Functor&& f) {
  auto pattern = details::getDefaultForEachExecutionPattern();
  for_each_index(pattern, range, std::forward<Functor>(f));
}

template <typename Functor>
void for_each_index(
  ExecutionPattern pattern,
  IndexSpace2D range,
  Functor&& f) {
  if (pattern == ExecutionPattern::parallel) {
      // Disable OpenMP parallelization for old Intel compilers to avoid compiler errors
#ifndef __INTEL_COMPILER
    #pragma omp parallel for collapse(2)
#endif
    for (atlas::idx_t i = range.i.start; i < range.i.end; ++i) {
      for (atlas::idx_t j = range.j.start; j < range.j.end; ++j) {
        f(i, j);
      }
    }
  } else if (pattern == ExecutionPattern::serial) {
    for (atlas::idx_t i = range.i.start; i < range.i.end; ++i) {
      for (atlas::idx_t j = range.j.start; j < range.j.end; ++j) {
        f(i, j);
      }
    }
  } else {
    throw eckit::BadParameter("Unknown execution pattern.", Here());
  }
}

template <typename Functor>
void for_each_index(
  IndexSpace2D range,
  Functor&& f) {
  auto pattern = details::getDefaultForEachExecutionPattern();
  for_each_index(pattern, range, std::forward<Functor>(f));
}

template <typename Functor>
void for_each_index(
  ExecutionPattern pattern,
  IndexSpace3D range,
  Functor&& f) {
  if (pattern == ExecutionPattern::parallel) {
      // Disable OpenMP parallelization for old Intel compilers to avoid compiler errors
#ifndef __INTEL_COMPILER
    #pragma omp parallel for collapse(3)
#endif
    for (atlas::idx_t i = range.i.start; i < range.i.end; ++i) {
      for (atlas::idx_t j = range.j.start; j < range.j.end; ++j) {
        for (atlas::idx_t k = range.k.start; k < range.k.end; ++k) {
          f(i, j, k);
        }
      }
    }
  } else if (pattern == ExecutionPattern::serial) {
    for (atlas::idx_t i = range.i.start; i < range.i.end; ++i) {
      for (atlas::idx_t j = range.j.start; j < range.j.end; ++j) {
        for (atlas::idx_t k = range.k.start; k < range.k.end; ++k) {
          f(i, j, k);
        }
      }
    }
  } else {
    throw eckit::BadParameter("Unknown execution pattern.", Here());
  }
}

template <typename Functor>
void for_each_index(
  IndexSpace3D range,
  Functor&& f) {
  auto pattern = details::getDefaultForEachExecutionPattern();
  for_each_index(pattern, range, std::forward<Functor>(f));
}

// The functor provided to this overload should have the following type:
//   (atlas::idx_t i, atlas::idx_t n, atlas::idx_t m) -> void
// where:
//   - i is the "flattened" index of the spectral coefficients. In Atlas recall that
//     spectral fields are stored in 2D arrays of "real" values where the first axis
//     is an ordering over the spectral coefficients and the second axis is the levels.
//     Furthermore, because the spectral coefficients are complex numbers, the real and imaginary
//     parts of each coefficient are stored in consecutive indices. The flattened index i is the
//     linear index into the parts of the coefficients. I.e. suppose `f` is a spectral field,
//     and i is even, then f(i, j) = Re(a^m_n) and f(i+1, j) = Im(a^m_n), where a^m_n is the
//     i/2th spectral coefficient on level j.
//   - n, the total wavenumber (degree) of the i/2th spectral coefficient.
//   - m, the zonal wavenumber (order)  of the i/2th spectral coefficient.
//
// This overload iterates over all spectral coefficients defined by the provided
// IndexSpaceSpectral1D, calling the functor for both the real and imaginary parts
// of each coefficient, passing the corresponding flattened index, n, and m. If you
// wish to restrict the iteration space, one can add appropriate if-conditions within
// the functor.
template <
  typename Functor,
  std::enable_if_t<
    std::is_invocable_v<
      Functor, atlas::idx_t, atlas::idx_t, atlas::idx_t>
    , int
  > = 0
>
void for_each_index(
  ExecutionPattern pattern,
  IndexSpaceSpectral1D range,
  Functor&& f) {
  if (pattern == ExecutionPattern::parallel) {
    details::ComputeSpectralCoefficientIndex compute_index(
      range.zonal_wavenumbers(), range.truncation());

    const int nb_zonal_wavenumbers{static_cast<int>(range.zonal_wavenumbers().size())};
      // Disable OpenMP parallelization for old Intel compilers to avoid compiler errors
#ifndef __INTEL_COMPILER
    #pragma omp parallel for
#endif
    for (int jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
      const int m = range.zonal_wavenumbers(jm);
      for (atlas::idx_t n = m; n <= range.truncation(); ++n) {
        auto[real_index, imag_index] = compute_index(n, jm, m);
        f(real_index, n, m);
        f(imag_index, n, m);
      }
    }
  } else if (pattern == ExecutionPattern::serial) {
    atlas::idx_t index = 0;
    const int nb_zonal_wavenumbers{static_cast<int>(range.zonal_wavenumbers().size())};
    for (int jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
        const int m = range.zonal_wavenumbers(jm);
        for (atlas::idx_t n = m; n <= range.truncation(); ++n) {
          f(index, n, m);
          f(index + 1, n, m);
          index += 2;
        }
    }
  } else {
    throw eckit::BadParameter("Unknown execution pattern.", Here());
  }
}

// The functor provided to this overload should have the following type:
//   (atlas::idx_t i_real, atlas::idx_t i_imag atlas::idx_t n, atlas::idx_t m) -> void
// where:
//   - i_real, the index of the real part of the i_real/2th spectral coefficient.
//   - i_imag, the index of the imaginary part of the i_real/2th spectral coefficient.
//             Note, i_imag == i_real + 1.
//   - n, the total wavenumber (degree) of the i_real/2th spectral coefficient.
//   - m, the zonal wavenumber (order)  of the i_real/2th spectral coefficient.
//
// This overload iterates over all spectral coefficients defined by the provided
// IndexSpaceSpectral1D, calling the functor for both the real and imaginary parts
// of each coefficient, passing the corresponding flattened i_real, i_imag, n, and m.
// If you wish to restrict the iteration space, one can add appropriate if-conditions
// within the functor.
template <
  typename Functor,
  std::enable_if_t<
    std::is_invocable_v<
      Functor, atlas::idx_t, atlas::idx_t, atlas::idx_t, atlas::idx_t>
    , int
  > = 0
>
void for_each_index(
  ExecutionPattern pattern,
  IndexSpaceSpectral1D range,
  Functor&& f) {
  if (pattern == ExecutionPattern::parallel) {
    details::ComputeSpectralCoefficientIndex compute_index(
      range.zonal_wavenumbers(), range.truncation());

    const int nb_zonal_wavenumbers{static_cast<int>(range.zonal_wavenumbers().size())};
      // Disable OpenMP parallelization for old Intel compilers to avoid compiler errors
#ifndef __INTEL_COMPILER
    #pragma omp parallel for
#endif
    for (int jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
      const int m = range.zonal_wavenumbers(jm);
      for (atlas::idx_t n = m; n <= range.truncation(); ++n) {
        auto[real_index, imag_index] = compute_index(n, jm, m);
        f(real_index, imag_index, n, m);
      }
    }
  } else if (pattern == ExecutionPattern::serial) {
    atlas::idx_t index = 0;
    const int nb_zonal_wavenumbers{static_cast<int>(range.zonal_wavenumbers().size())};
    for (int jm = 0; jm < nb_zonal_wavenumbers; ++jm) {
      const int m = range.zonal_wavenumbers(jm);
      for (atlas::idx_t n = m; n <= range.truncation(); ++n) {
        f(index, index + 1, n, m);
        index += 2;
      }
    }
  } else {
    throw eckit::BadParameter("Unknown execution pattern.", Here());
  }
}

// The functor provided to this overload should have one of the following types:
// (atlas::idx_t i, atlas::idx_t n, atlas::idx_t m) -> void
// (atlas::idx_t i_real, atlas::idx_t i_imag, atlas::idx_t n, atlas::idx_t m) -> void
// See the above overloads for details on the parameters.
template <typename Functor>
void for_each_index(
  IndexSpaceSpectral1D range,
  Functor&& f) {
  auto pattern = details::getDefaultForEachExecutionPattern();
  for_each_index(pattern, range, std::forward<Functor>(f));
}

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
    // When owned points are not contiguous, we must examine ghost field
    const auto ghost = atlas::array::make_view<int, 1>(firstField.functionspace().ghost());
    const atlas::idx_t iMax = firstField.shape(0);
    const atlas::idx_t jMax = firstField.shape(1);
    if (pattern == ExecutionPattern::parallel) {
      // Disable OpenMP parallelization for old Intel compilers to avoid compiler errors
#ifndef __INTEL_COMPILER
      #pragma omp parallel for
#endif
      for (atlas::idx_t i = 0; i < iMax; ++i) {
        if (ghost(i) == 0) {
          for (atlas::idx_t j = 0; j < jMax; ++j) {
            std::apply([&](auto&&... views) { f(views(i, j)...); }, viewsTuple);
          }
        }
      }
    } else if (pattern == ExecutionPattern::serial) {
      for (atlas::idx_t i = 0; i < iMax; ++i) {
        if (ghost(i) == 0) {
          for (atlas::idx_t j = 0; j < jMax; ++j) {
            std::apply([&](auto&&... views) { f(views(i, j)...); }, viewsTuple);
          }
        }
      }
    } else {
      throw eckit::BadParameter("Unknown execution pattern.", Here());
    }
  } else {
    const bool include_halo = range == IndexRange::include_halo;
    const auto range_pair = details::getContiguousHorizontalIndexSpace(firstField, include_halo);
    const atlas::idx_t iStart = range_pair.first;
    const atlas::idx_t iEnd = range_pair.second;
    const atlas::idx_t jMax = firstField.shape(1);
    if (pattern == ExecutionPattern::parallel) {
      // Disable OpenMP parallelization for old Intel compilers to avoid compiler errors
#ifndef __INTEL_COMPILER
      #pragma omp parallel for collapse(2)
#endif
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
    // When owned points are not contiguous, we must examine ghost field
    const auto ghost = atlas::array::make_view<int, 1>(firstField.functionspace().ghost());
    const atlas::idx_t iMax = firstField.shape(0);
    if (pattern == ExecutionPattern::parallel) {
      // Disable OpenMP parallelization for old Intel compilers to avoid compiler errors
#ifndef __INTEL_COMPILER
      #pragma omp parallel for
#endif
      for (atlas::idx_t i = 0; i < iMax; ++i) {
        if (ghost(i) == 0) {
          std::apply([&](auto&&... views) {
            f(views.slice(i, atlas::array::Range::all())...);
          }, viewsTuple);
        }
      }
    } else if (pattern == ExecutionPattern::serial) {
      for (atlas::idx_t i = 0; i < iMax; ++i) {
        if (ghost(i) == 0) {
          std::apply([&](auto&&... views) {
            f(views.slice(i, atlas::array::Range::all())...);
          }, viewsTuple);
        }
      }
    } else {
      throw eckit::BadParameter("Unknown execution pattern.", Here());
    }
  } else {
    const bool include_halo = range == IndexRange::include_halo;
    const auto range_pair = details::getContiguousHorizontalIndexSpace(firstField, include_halo);
    const atlas::idx_t iStart = range_pair.first;
    const atlas::idx_t iEnd = range_pair.second;
    if (pattern == ExecutionPattern::parallel) {
      // Disable OpenMP parallelization for old Intel compilers to avoid compiler errors
#ifndef __INTEL_COMPILER
      #pragma omp parallel for
#endif
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

// perThreadStorage: a convience routine for allocating a buffer for use as per thread
// storage within a parallel region e.g. for_each_value, for_each_column, for_each_index.
// Let `t` be the number of threads, and `(s1, s2,..., sn)` be the user provided per
// thread array shape, then the returned array has shape `(t, s1, s2,..., sn)`. Within the
// parallel region, the buffer view can be sliced to provides a specific threads storage
// using `atlas_omp_get_thread_num()`. For example:
// ```
// auto buffer = perThreadStorage<double>(k);
// auto bufferView = array::make_view<double, 2>(buffer);
// for_each_index(
//   index_space_1d,
//   [=](idx_t i) mutable {
//     ...
//     auto threadChunkView = bufferView.slice(
//       atlas_omp_get_thread_num(), atlas::array::Range::all());
//     ...
//   });
// ```
template<typename T, typename... S>
atlas::array::ArrayT<T> perThreadStorage(
      util::ExecutionPattern pattern, S... per_thread_shape)
{
  if (pattern == util::ExecutionPattern::parallel) {
      return atlas::array::ArrayT<T>(
          static_cast<atlas::idx_t>(atlas_omp_get_max_threads()), per_thread_shape...);
  } else {
      return atlas::array::ArrayT<T>(1,  per_thread_shape...);
  }
};

// Overload for default execution pattern
template<typename T, typename... S>
atlas::array::ArrayT<T> perThreadStorage(S... per_thread_shape)
{
  return perThreadStorage<T>(util::details::getDefaultForEachExecutionPattern(),
    per_thread_shape...);
};

}  // namespace util
