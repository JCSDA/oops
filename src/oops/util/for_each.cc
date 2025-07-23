/*
 * (C) Copyright 2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#include "oops/util/for_each.h"

namespace util {
namespace details {

ExecutionPattern getDefaultForEachExecutionPattern() { return ExecutionPattern::parallel; }

IndexRange getDefaultForEachIndexRange() { return IndexRange::exclude_halo; }

}  // namespace details

IndexSpace1D make_index_space_1d(IndexRange range, const atlas::Field& field) {
  if (range == IndexRange::exclude_halo && !details::hasContiguousOwnedPoints(field)) {
    throw eckit::NotImplemented(
        "Excluding halo values is not implemented for non-contiguous halos.", Here());
  } else {
    const bool include_halo = range == IndexRange::include_halo;
    const auto[iStart, iEnd] = details::getContiguousHorizontalIndexSpace(field, include_halo);
    return IndexSpace1D{ iStart, iEnd };
  }
}

IndexSpace2D make_index_space_2d(IndexRange range, const atlas::Field& field) {
  return IndexSpace2D{
    make_index_space_1d(range, field),
    IndexSpace1D{0, field.shape(1)}
  };
}

IndexSpace3D make_index_space_3d(IndexRange range, const atlas::Field& field) {
  return IndexSpace3D{
    make_index_space_1d(range, field),
    IndexSpace1D{0, field.shape(1)},
    IndexSpace1D{0, field.shape(2)}
  };
}

IndexSpaceSpectral1D make_index_space_spectral_1d(const atlas::Field& field) {
  if (auto sp = atlas::functionspace::Spectral(field.functionspace())) {
    return IndexSpaceSpectral1D{ sp };
  } else {
    throw eckit::BadParameter("Not a Spectral FunctionSpace", Here());
  }
}

}  // namespace util
