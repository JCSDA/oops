// (C) Crown Copyright 2026 Met Office
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#include "oops/generic/InterpolationCache.h"

#include <algorithm>
#include <cstddef>
#include <iterator>

#include "atlas/array/MakeView.h"
#include "atlas/field/Field.h"
#include "atlas/functionspace/PointCloud.h"
#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace oops {

InterpolationCache::InterpolationCache(
    const atlas::FunctionSpace& sourceFunctionSpace,
    const eckit::LocalConfiguration& interpMethod,
    const std::vector<double>& targetLats,
    const std::vector<double>& targetLons)
    : sourceFunctionSpace_{sourceFunctionSpace},
      interpMethod_{[&] {
        auto subConf = interpMethod;
        if (subConf.getBool("allow_halo_exchange", false)) {
          throw eckit::BadValue(
              "InterpolationCache: allow_halo_exchange must be false.", Here());
        }
        subConf.set("allow_halo_exchange", false);
        return subConf;
      }()},
      targetLonLats_{[&] {
        if (targetLats.size() != targetLons.size()) {
          throw eckit::BadValue(
              "Target lat and lon vectors must be the same size.", Here());
        }
        auto lonLats = std::vector<atlas::PointLonLat>{};
        lonLats.reserve(targetLats.size());
        std::transform(targetLats.begin(), targetLats.end(), targetLons.begin(),
                       std::back_inserter(lonLats),
                       [](const double lat, const double lon) {
                         auto point = atlas::PointLonLat{lon, lat};
                         point.normalise();
                         return point;
                       });
        return lonLats;
      }()} {
  Log::trace() << classname() + "::InterpolationCache start" << std::endl;
  Log::trace() << classname() + "::InterpolationCache done" << std::endl;
}

const atlas::Interpolation& InterpolationCache::getInterp(
    const std::vector<bool>& mask) const {
  Log::trace() << classname() + "::getInterp start" << std::endl;
  util::Timer timer(classname(), "getInterp");

  if (mask.size() != targetLonLats_.size()) {
    throw eckit::BadValue(
        "Mask size must match the number of target lon-lat points.", Here());
  }

  // Find or insert interpolation object in map.
  auto& interp = interpMap_[mask];

  if (!interp) {
    // Make a lon-lat field.
    const atlas::idx_t fieldSize = std::count(mask.cbegin(), mask.cend(), true);
    auto lonLatField =
        atlas::Field("lonlat", atlas::array::make_datatype<double>(),
                     atlas::array::make_shape(fieldSize, 2));

    auto lonLatView = atlas::array::make_view<double, 2>(lonLatField);

    // Copy lonlats with mask == true.
    atlas::idx_t maskedIdx = 0;
    for (std::size_t unmaskedIdx = 0; unmaskedIdx < mask.size();
         unmaskedIdx++) {
      if (mask[unmaskedIdx] == true) {
        lonLatView(maskedIdx, 0) = targetLonLats_[unmaskedIdx].lon();
        lonLatView(maskedIdx, 1) = targetLonLats_[unmaskedIdx].lat();
        maskedIdx++;
      }
    }

    // Create an interpolation object.
    const auto targetFunctionSpace =
        atlas::functionspace::PointCloud(lonLatField);
    interp = atlas::Interpolation(interpMethod_, sourceFunctionSpace_,
                                  targetFunctionSpace);
  }

  Log::trace() << classname() + "::getInterp done" << std::endl;
  return interp;
}

}  // namespace oops
