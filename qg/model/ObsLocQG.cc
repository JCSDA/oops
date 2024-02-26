/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/ObsLocQG.h"

#include <memory>

#include "atlas/array.h"
#include "eckit/config/Configuration.h"
#include "eckit/geometry/Point2.h"
#include "eckit/geometry/Point3.h"
#include "eckit/geometry/Sphere.h"

#include "model/GeometryQGIterator.h"
#include "model/LocationsQG.h"
#include "model/ObsSpaceQG.h"
#include "model/QgTraits.h"

using atlas::array::make_view;

namespace qg {

static oops::ObsLocalizationMaker<QgTraits, QgObsTraits, ObsLocQG> makerObsLoc_("Heaviside");

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

  for (int jj = 0; jj < locs->size(); ++jj) {
    eckit::geometry::Point2 obsPoint(lonlat(jj, 0), lonlat(jj, 1));
    double localDist = eckit::geometry::Sphere::distance(6.371e6, refPoint2, obsPoint);
    if (localDist > lengthscale_) {
      local.setToMissing(jj);
    } else {
      local *= 1.0;
    }
  }
}

// -----------------------------------------------------------------------------

void ObsLocQG::print(std::ostream & os) const {
  os << "Observation space localization: Heaviside with lengthscale = " << lengthscale_;
}

// -----------------------------------------------------------------------------

}  // namespace qg
