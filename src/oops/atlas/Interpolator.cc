/*
 * (C) Copyright 2023-2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/atlas/Interpolator.h"

#include <string>

#include "atlas/array.h"
#include "atlas/field.h"

#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"

namespace oops {
namespace atlasbase {

// -----------------------------------------------------------------------------

Interpolator::Interpolator(const eckit::Configuration & /*conf*/,
                           const GeometryData & /*source_geometry*/,
                           const std::vector<double> & target_lats,
                           const std::vector<double> & target_lons) {
  ASSERT(target_lats.size() == target_lons.size());
  nb_targets_ = target_lats.size();
}

// -----------------------------------------------------------------------------

void Interpolator::apply(const Variables& vars,
                         const atlas::FieldSet& fields,
                         std::vector<double>& buffer) const {
  std::vector<bool> mask(nb_targets_, true);
  apply(vars, fields, mask, buffer);
}

// -----------------------------------------------------------------------------

void Interpolator::applyAD(const Variables& vars,
                           atlas::FieldSet& fields,
                           const std::vector<double>& buffer) const {
  std::vector<bool> mask(nb_targets_, true);
  applyAD(vars, fields, mask, buffer);
}

// -----------------------------------------------------------------------------

void Interpolator::preprocess(atlas::FieldSet & fields) {
  fields.haloExchange();
}

// -----------------------------------------------------------------------------

void Interpolator::preprocessAD(atlas::FieldSet & fields) {
  fields.adjointHaloExchange();
  fields.set_dirty();
}

// -----------------------------------------------------------------------------

// Unscramble MPI buffer into the model's FieldSet representation
void Interpolator::bufferToFieldSet(const Variables & vars,
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
    for (size_t ji = 0; ji < buffer_chunk_size; ++ji) {
      for (size_t jlev = 0; jlev < num_levels; ++jlev, ++current) {
        const size_t index = buffer_indices[ji];
        ASSERT(static_cast<size_t>(std::distance(buffer_start, current)) < buffer_size);
        view(index, jlev) = *current;
      }
    }
  }
}

// -----------------------------------------------------------------------------

// (Adjoint of) Unscramble MPI buffer into the model's FieldSet representation
void Interpolator::bufferToFieldSetAD(const Variables & vars,
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
    for (size_t ji = 0; ji < buffer_chunk_size; ++ji) {
      for (size_t jlev = 0; jlev < num_levels; ++jlev, ++current) {
        const size_t index = buffer_indices[ji];
        ASSERT(static_cast<size_t>(std::distance(buffer_start, current)) < buffer_size);
        *current += view(index, jlev);
      }
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace atlasbase
}  // namespace oops
