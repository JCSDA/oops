/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "lorenz95/ObsLocGC99.h"

#include <algorithm>
#include <cmath>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/geometry/Point2.h"

#include "lorenz95/Iterator.h"
#include "lorenz95/L95Traits.h"
#include "lorenz95/ObsTable.h"
#include "lorenz95/ObsVec1D.h"

#include "oops/generic/gc99.h"
#include "oops/interface/ObsLocalization.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {

static oops::ObsLocalizationMaker<L95Traits, L95ObsTraits,
              oops::ObsLocalization<L95Traits, L95ObsTraits, ObsLocGC99>> makerGC_("Gaspari-Cohn");

// -----------------------------------------------------------------------------

ObsLocGC99::ObsLocGC99(const eckit::Configuration & config, const ObsTable & obsdb)
  : rscale_(config.getDouble("lengthscale")), obsdb_(obsdb)
{
}

// -----------------------------------------------------------------------------

void ObsLocGC99::computeLocalization(const Iterator & iterator, ObsData1D<int> & outside,
                                     ObsVec1D & result) const {
  std::vector<double> locations = obsdb_.locations();
  eckit::geometry::Point2 center = *iterator;
  for (unsigned int ii=0; ii < obsdb_.nobs(); ++ii) {
    double curdist = std::abs(center[0] - locations[ii]);
    curdist = std::min(curdist, 1.-curdist);
    if (curdist >= rscale_) {
      outside[ii] = 1;
    } else {
      outside[ii] = 0;
      result[ii] = oops::gc99(curdist/rscale_);
    }
  }
}

// -----------------------------------------------------------------------------

void ObsLocGC99::print(std::ostream & os) const {
  os << "Gaspari-Cohn localization with lengthscale=" << rscale_;
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
