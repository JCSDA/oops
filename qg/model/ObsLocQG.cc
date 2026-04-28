/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/ObsLocQG.h"

#include <memory>
#include <vector>

#include "atlas/array.h"
#include "eckit/config/Configuration.h"
#include "eckit/geometry/Point2.h"
#include "eckit/geometry/Point3.h"
#include "eckit/geometry/Sphere.h"

#include "model/GeometryQGIterator.h"
#include "model/LocationsQG.h"
#include "model/ObsSpaceQG.h"
#include "model/ObsVecQG.h"
#include "oops/generic/gc99.h"


using atlas::array::make_view;

namespace qg {

// -----------------------------------------------------------------------------

ObsLocQG::ObsLocQG(const eckit::Configuration & conf, const ObsSpaceQG & obsdb)
  : lengthscale_(conf.getDouble("lengthscale")), obsdb_(obsdb)
{
}

// -----------------------------------------------------------------------------

void ObsLocQG::computeLocalization(const GeometryQGIterator & p,
                                   ObsVecQG & local) const {
  std::unique_ptr<LocationsQG> locs = obsdb_.locations();
  atlas::Field field_lonlat = locs->lonlat();
  auto lonlat = make_view<double, 2>(field_lonlat);
  eckit::geometry::Point3 refPoint = *p;
  eckit::geometry::Point2 refPoint2(refPoint[0], refPoint[1]);

  // get the number of variables per location.
  // This feels a bit hacky, but it is not clear how to get nvar otherwise
  const size_t nvar = local.size() / locs->size();

  // Calculate the localization weights for all locations and variables
  ObsVecQG weights(local);
  weights.ones();
  std::vector<double> vec;
  weights.serialize(vec);
  for (int jj = 0; jj < locs->size(); ++jj) {
    // calculate the localization weight for this location
    const eckit::geometry::Point2 obsPoint(lonlat(jj, 0), lonlat(jj, 1));
    const double localDist = eckit::geometry::Sphere::distance(6.371e6, refPoint2, obsPoint);
    const double weight = oops::gc99(localDist / lengthscale_);
    for (int kk = 0; kk < nvar; ++kk) {
      const size_t idx = jj * nvar + kk;
      if (localDist > lengthscale_) {
        vec[idx] = 0.0;
      } else {
        vec[idx] *= weight;
      }
    }
  }

  // deserialize the weights back into the ObsVecQG
  size_t s = 0;
  weights.deserialize(vec, s);

  // set to missing value if the weight is zero
  for (size_t jj = 0; jj < locs->size(); ++jj) {
    if (vec[jj * nvar] <= 0.0) {
      weights.setToMissing(jj);
    }
  }

  // finally, multiply the input local by the weights
  local *= weights;
}

// -----------------------------------------------------------------------------

double ObsLocQG::computeLocalization(const eckit::geometry::Point3 & point1,
                                     const eckit::geometry::Point3 & point2) const {
  eckit::geometry::Point2 point1_2(point1[0], point1[1]);
  eckit::geometry::Point2 point2_2(point2[0], point2[1]);
  double localDist = eckit::geometry::Sphere::distance(6.371e6, point1_2, point2_2);
  if (localDist > lengthscale_) {
    return 0.0;
  } else {
    return oops::gc99(localDist / lengthscale_);
  }
}

// -----------------------------------------------------------------------------

void ObsLocQG::print(std::ostream & os) const {
  os << "Observation space localization: GC99 with lengthscale = " << lengthscale_;
}

// -----------------------------------------------------------------------------

}  // namespace qg
