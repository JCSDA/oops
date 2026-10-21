// (C) Crown Copyright 2026 Met Office
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "atlas/functionspace/FunctionSpace.h"
#include "atlas/interpolation/Interpolation.h"
#include "atlas/util/Point.h"
#include "eckit/config/LocalConfiguration.h"

namespace oops {

/// \brief Lazily builds and caches atlas::Interpolation objects, keyed by
/// location mask.
///
/// Each distinct mask over the target lon-lats requires a different (and
/// expensive to construct) atlas::Interpolation object. This class defers
/// construction until a mask is first requested, and reuses the result for any
/// subsequent request with the same mask.
class InterpolationCache {
 public:
  /// \brief Class name string.
  static std::string classname() { return "oops::InterpolationCache"; }

  /// \brief Construct a cache of interpolations from \p sourceFunctionSpace to
  /// the target points given by \p targetLats/\p targetLons, using \p
  /// interpMethod.
  ///
  /// \p interpMethod must not set allow_halo_exchange to true, since the
  /// source fields are assumed to have had halo exchanges performed beforehand,
  /// and adjoints thereof afterwards.
  InterpolationCache(const atlas::FunctionSpace& sourceFunctionSpace,
                     const eckit::LocalConfiguration& interpMethod,
                     const std::vector<double>& targetLats,
                     const std::vector<double>& targetLons);

  /// \brief Get the interpolation object for \p mask, creating and caching it
  /// if necessary.
  const atlas::Interpolation& getInterp(const std::vector<bool>& mask) const;

 private:
  // FunctionSpace from Geometry.
  atlas::FunctionSpace sourceFunctionSpace_;

  // Atlas Interpolation method config.
  eckit::LocalConfiguration interpMethod_;

  // Vector of lon-lats from constructor.
  std::vector<atlas::PointLonLat> targetLonLats_{};

  // Map of interpolation objects, keyed by mask.
  // Mutable keyword required as we need to construct an interpolation object if
  // it doesn't already exist.
  mutable std::unordered_map<std::vector<bool>, atlas::Interpolation>
      interpMap_{};
};

}  // namespace oops
