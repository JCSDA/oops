/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "lorenz95/ObsLocBoxCar.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/geometry/Point3.h"

#include "lorenz95/Iterator.h"
#include "lorenz95/L95Traits.h"
#include "lorenz95/ObsTable.h"
#include "lorenz95/ObsVec1D.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {

static oops::ObsLocalizationMaker<L95Traits, L95ObsTraits, ObsLocBoxCar> makerGC_("Box Car");

// -----------------------------------------------------------------------------

ObsLocBoxCar::ObsLocBoxCar(const eckit::Configuration & config, const ObsTable & obsdb)
  : rscale_(config.getDouble("lengthscale")), obsdb_(obsdb)
{}

// -----------------------------------------------------------------------------

void ObsLocBoxCar::computeLocalization(const Iterator & iterator, ObsVec1D & locfactor) const {
  std::vector<double> locations = obsdb_.locations();
  eckit::geometry::Point3 center = *iterator;
  for (unsigned int ii=0; ii < obsdb_.nobs(); ++ii) {
    double curdist = std::abs(center[0] - locations[ii]);
    curdist = std::min(curdist, 1.-curdist);
    if (curdist >= rscale_) {
      locfactor[ii] = locfactor.missing();
    } else {
      locfactor[ii] *= 1.0;
    }
  }
}

// -----------------------------------------------------------------------------

void ObsLocBoxCar::print(std::ostream & os) const {
  os << "Box Car localization with lengthscale=" << rscale_;
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
