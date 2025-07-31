/*
 * (C) Copyright 2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#include "oops/util/reduction.h"

#include <new>
#include <vector>

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/FunctionSpace.h"
#include "oops/util/for_each.h"
#include "oops/util/missingValues.h"

namespace util {

namespace details {

ExecutionPattern getDefaultDotProductExecutionPattern() { return ExecutionPattern::parallel; }

#if defined(__cpp_lib_hardware_interference_size)
constexpr int cache_line_size = std::hardware_destructive_interference_size;
#else
constexpr int cache_line_size = 64;
#endif

struct CacheLineDouble {
  alignas(cache_line_size) double value = 0.0;
};

std::vector<CacheLineDouble> get_reduction_buffer(const ExecutionPattern pattern) {
  if (pattern == util::ExecutionPattern::parallel) {
    return std::vector<CacheLineDouble>(atlas_omp_get_max_threads());
  } else {
    return std::vector<CacheLineDouble>(1);
  }
}

inline double product_unless_missing(const double x1, const double x2) {
  static const double missing = missingValue<double>();
  if (x1 != missing && x2 != missing) {
    return x1 * x2;
  } else {
    return 0.0;
  }
}

}  // namespace details

double dot_product_spectral(const ExecutionPattern pattern,
    const atlas::Field & field1, const atlas::Field & field2) {
  const auto view1 = atlas::array::make_view<double, 2>(field1);
  const auto view2 = atlas::array::make_view<double, 2>(field2);

  // Compute dot_product per-thread, write result into buffer indexed by thread number
  auto buffer = details::get_reduction_buffer(pattern);
  for_each_index(
    pattern,
    make_index_space_spectral_1d(field1),
    [&](atlas::idx_t real_i, atlas::idx_t imag_i, atlas::idx_t n, atlas::idx_t m) {
      const double real_fact = (m == 0) ? 1.0 : 2.0;
      const double imag_fact = (m == 0) ? 0.0 : 2.0;
      for (atlas::idx_t j = 0; j < field1.shape(1); ++j) {
        buffer[atlas_omp_get_thread_num()].value +=
            real_fact * details::product_unless_missing(view1(real_i, j), view2(real_i, j));
        buffer[atlas_omp_get_thread_num()].value +=
            imag_fact * details::product_unless_missing(view1(imag_i, j), view2(imag_i, j));
      }
    });

  // Reduce over the buffer (in serial: it's too short to be worth multithreading)
  double sum = 0.0;
  for (size_t i = 0; i < buffer.size(); ++i) {
    sum += buffer[i].value;
  }

  return sum;
}

double dot_product_nodal(const ExecutionPattern pattern,
    const atlas::Field & field1, const atlas::Field & field2) {
  // Compute dot_product per-thread, write result into buffer indexed by thread number
  auto buffer = details::get_reduction_buffer(pattern);
  for_each_value(
    pattern,
    util::IndexRange::exclude_halo,
    [&](const double val1, const double val2) {
      buffer[atlas_omp_get_thread_num()].value += details::product_unless_missing(val1, val2);
    },
    field1, field2);

  // Reduce over the buffer (in serial: it's too short to be worth multithreading)
  double sum = 0.0;
  for (size_t i = 0; i < buffer.size(); ++i) {
    sum += buffer[i].value;
  }

  return sum;
}

double dot_product_on_task(const ExecutionPattern pattern,
    const atlas::Field & field1, const atlas::Field & field2) {
  ASSERT(field1.rank() == 2);
  ASSERT(field2.rank() == 2);
  ASSERT(field1.shape(0) == field2.shape(0));
  ASSERT(field1.shape(1) == field2.shape(1));
  ASSERT(field1.functionspace().type() == field2.functionspace().type());

  if (field1.functionspace().type() == "Spectral") {
    return dot_product_spectral(pattern, field1, field2);
  } else {
    return dot_product_nodal(pattern, field1, field2);
  }
}

double dot_product_on_task(const atlas::Field & field1, const atlas::Field & field2) {
  auto pattern = details::getDefaultDotProductExecutionPattern();
  return dot_product_on_task(pattern, field1, field2);
}

}  // namespace util
