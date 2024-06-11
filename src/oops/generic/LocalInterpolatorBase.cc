/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/generic/LocalInterpolatorBase.h"

#include <string>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"

namespace oops {

// Unscramble MPI buffer into the model's FieldSet representation
void LocalInterpolatorBase::bufferToFieldSet(const Variables & vars,
                                             const std::vector<size_t> & buffer_indices,
                                             const std::vector<double> & buffer,
                                             atlas::FieldSet & target) {
  const size_t buffer_chunk_size = buffer_indices.size();
  const size_t buffer_size = buffer.size();
  ASSERT(buffer_chunk_size > 0);
  ASSERT(buffer_size % buffer_chunk_size == 0);

  const auto buffer_start = buffer.begin();
  auto current = buffer.begin();

  for (size_t jf = 0; jf < vars.size(); ++jf) {
    const std::string & fname = vars[jf].name();
    atlas::Field & field = target.field(fname);

    atlas::array::ArrayView<double, 2> view = atlas::array::make_view<double, 2>(field);
    const size_t field_size = view.shape(0);
    const size_t num_levels = view.shape(1);
    ASSERT(buffer_chunk_size <= field_size);
    for (size_t jlev = 0; jlev < num_levels; ++jlev) {
      for (size_t ji = 0; ji < buffer_chunk_size; ++ji, ++current) {
        const size_t index = buffer_indices[ji];
        ASSERT(static_cast<size_t>(std::distance(buffer_start, current)) < buffer_size);
        view(index, jlev) = *current;
      }
    }
  }
}

// -----------------------------------------------------------------------------

// (Adjoint of) Unscramble MPI buffer into the model's FieldSet representation
void LocalInterpolatorBase::bufferToFieldSetAD(const Variables & vars,
                                               const std::vector<size_t> & buffer_indices,
                                               std::vector<double> & buffer,
                                               const atlas::FieldSet & target) {
  const size_t buffer_chunk_size = buffer_indices.size();
  const size_t buffer_size = buffer.size();
  ASSERT(buffer_chunk_size > 0);
  ASSERT(buffer_size % buffer_chunk_size == 0);

  const auto buffer_start = buffer.begin();
  auto current = buffer.begin();

  for (size_t jf = 0; jf < vars.size(); ++jf) {
    const std::string & fname = vars[jf].name();
    atlas::Field & field = target.field(fname);  // const in principle, but intel can't compile that

    const atlas::array::ArrayView<double, 2> view = atlas::array::make_view<double, 2>(field);
    const size_t field_size = view.shape(0);
    const size_t num_levels = view.shape(1);
    ASSERT(buffer_chunk_size <= field_size);
    for (size_t jlev = 0; jlev < num_levels; ++jlev) {
      for (size_t ji = 0; ji < buffer_chunk_size; ++ji, ++current) {
        const size_t index = buffer_indices[ji];
        ASSERT(static_cast<size_t>(std::distance(buffer_start, current)) < buffer_size);
        *current += view(index, jlev);
      }
    }
  }
}

}  // namespace oops
