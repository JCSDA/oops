/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "lorenz95/ObsLocGC99.h"

#include <vector>

#include "eckit/config/Configuration.h"

#include "lorenz95/L95Traits.h"
#include "lorenz95/ObsTableView.h"
#include "lorenz95/ObsVec1D.h"

#include "oops/generic/gc99.h"
#include "oops/interface/ObsLocalization.h"
#include "oops/util/missingValues.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {

static oops::ObsLocalizationMaker<L95Traits, L95ObsTraits,
              oops::ObsLocalization<L95Traits, L95ObsTraits, ObsLocGC99>> makerGC_("Gaspari-Cohn");

// -----------------------------------------------------------------------------

ObsLocGC99::ObsLocGC99(const eckit::Configuration & config, const ObsTableView & obsdb)
  : rscale_(config.getDouble("lengthscale")), obsdb_(obsdb)
{
}

// -----------------------------------------------------------------------------

void ObsLocGC99::computeLocalization(const Iterator &, ObsVec1D & result) const {
  const std::vector<double> & obsdist = obsdb_.obsdist();
  double missing = util::missingValue(missing);
  for (unsigned int ii=0; ii < obsdb_.nobs(); ++ii) {
    result[ii] = oops::gc99(obsdist[ii]/rscale_);
  }
}

// -----------------------------------------------------------------------------

void ObsLocGC99::print(std::ostream & os) const {
  os << "Gaspari-Cohn localization with lengthscale=" << rscale_;
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
