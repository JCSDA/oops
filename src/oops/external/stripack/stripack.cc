/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/external/stripack/stripack.h"

#include <algorithm>
#include <array>
#include <vector>

#include "atlas/util/Geometry.h"
#include "atlas/util/Point.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/FloatCompare.h"
#include "oops/util/Random.h"

namespace stripack {

Triangulation::Triangulation(const std::vector<double> & lats, const std::vector<double> & lons)
  : num_nodes_(lats.size()),
    xs_(num_nodes_),
    ys_(num_nodes_),
    zs_(num_nodes_),
    random_permutation_(num_nodes_),
    inverse_random_permutation_(num_nodes_),
    list_(6*(num_nodes_-2)),
    lptr_(6*(num_nodes_-2)),
    lend_(num_nodes_),
    near_(num_nodes_, 0),
    next_(num_nodes_, 0),
    dist_(num_nodes_, 0.0)
{
  ASSERT(lats.size() == lons.size());

  const atlas::Geometry unitsphere(1.0);

  // The source latlons come from a model grid so may be highly structured. But triangulating by
  // traversing through structured points can lead to a very non-uniform intermediate triangulation
  // with pathological configurations. So we randomize the source coordinates before triangulation.
  //
  // First, construct a random permutation and its inverse:
  std::vector<size_t> indices(num_nodes_);
  std::iota(indices.begin(), indices.end(), 0);
  const unsigned int seed = 7452;   // arbitrary, but reproducible
  util::shuffle(begin(indices), end(indices), seed);
  for (size_t i = 0; i < num_nodes_; ++i) {
    random_permutation_[i] = indices[i];
    inverse_random_permutation_[indices[i]] = i;
  }

  // Then permute source points while converting from latlon coords to xyz coords
  for (size_t i = 0; i < num_nodes_; ++i) {
    const size_t ip = inverse_random_permutation_[i];
    atlas::PointLonLat ptll(lons[ip], lats[ip]);
    ptll.normalise();
    atlas::Point3 pt3;
    unitsphere.lonlat2xyz(ptll, pt3);
    xs_[i] = pt3[0];
    ys_[i] = pt3[1];
    zs_[i] = pt3[2];
  }

  stripack_trmesh_f90(
      static_cast<int>(num_nodes_),
      xs_.data(), ys_.data(), zs_.data(),
      list_.data(), lptr_.data(), lend_.data(), lnew_,
      near_.data(), next_.data(), dist_.data());
}

bool Triangulation::containingTriangleAndBarycentricCoords(
      const std::array<double, 3> & coords,
      const int guess_index,
      std::array<int, 3> & vertexIndices,
      std::array<double, 3> & barycentricCoords) const {
  ASSERT(guess_index >= 0 && guess_index <= static_cast<int>(num_nodes_) - 1);
  ASSERT(oops::is_close_absolute(
      sqrt(coords[0]*coords[0] + coords[1]*coords[1] + coords[2]*coords[2]),
      1.0, 1.0e-9));

  // Randomize inputs
  const int randomized_guess_index = random_permutation_[guess_index];

  std::array<int, 3> tmp{};
  stripack_trfind_f90(
      randomized_guess_index,
      coords.data(),
      static_cast<int>(num_nodes_),
      xs_.data(), ys_.data(), zs_.data(),
      list_.data(), lptr_.data(), lend_.data(),
      barycentricCoords.data(),
      tmp.data());

  // If all indices are -1 (Fortran index 0), STRIPACK identified all points as coplanar. We can't
  // recover from this, so we ASSERT for this case.
  // If one (the last) index is -1 (Fortran index 0), STRIPACK identified the target point as being
  // outside the convex hull of the source points. We can recover from this by passing the invalid
  // index to the interpolator, which will in turn fill the target result with a missing value.
  constexpr int invalid = -1;
  const size_t numInvalid = std::count(tmp.cbegin(), tmp.cend(), invalid);
  ASSERT(numInvalid == 0 || numInvalid == 1);  // point is enclosed or is outside hull
  if (numInvalid == 1) {
    return false;
  }

  // Invert the randomization, so the returned index matches the coordinates used in construction
  for (size_t i = 0; i < 3; ++i) {
    vertexIndices[i] = inverse_random_permutation_[tmp[i]];
  }
  return true;
}

}  // namespace stripack
