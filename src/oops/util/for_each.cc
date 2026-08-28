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

ComputeSpectralCoefficientIndex::ComputeSpectralCoefficientIndex(
    atlas::array::LocalView<const int, 1> zonal_wavenumbers,
    int truncation)
  : prefix_sum(zonal_wavenumbers.size() + 1),
    prefix_sum_view(atlas::array::make_view<int, 1>(prefix_sum)) {
  prefix_sum_view.assign(0);
  for (int jm = 0; jm < zonal_wavenumbers.size(); ++jm) {
    prefix_sum_view(jm + 1) = prefix_sum_view(jm) + ((truncation + 1) - zonal_wavenumbers(jm));
  }
}

std::pair<atlas::idx_t, atlas::idx_t>
ComputeSpectralCoefficientIndex::operator()(int n, int jm, int m) {
  auto index = prefix_sum_view[jm] + (n - m);
  auto real_index = index * 2;
  auto imag_index = index * 2 + 1;
  return std::make_pair(real_index, imag_index);
}

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
