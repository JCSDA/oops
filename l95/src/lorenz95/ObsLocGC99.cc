/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "lorenz95/ObsLocGC99.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "lorenz95/ObsTableView.h"
#include "lorenz95/ObsVec1D.h"

#include "oops/generic/gc99.h"
#include "oops/interface/ObsLocalization.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {
static oops::ObsLocalizationMaker<L95Traits,
              oops::ObsLocalization<L95Traits, ObsLocGC99>> makerGC_("Gaspari-Cohn");

// -----------------------------------------------------------------------------

ObsLocGC99::ObsLocGC99(const eckit::Configuration & config, const ObsTableView & obsdb)
  : obsdb_(obsdb),
    rscale_(config.getDouble("lengthscale"))
{
}

// -----------------------------------------------------------------------------

ObsLocGC99::~ObsLocGC99() {}

// -----------------------------------------------------------------------------

void ObsLocGC99::multiply(ObsVec1D & dy) const {
  const std::vector<double> & obsdist = obsdb_.obsdist();
  double  gc;
  for (unsigned int ii=0; ii < dy.nobs(); ++ii) {
    gc = oops::gc99(obsdist[ii]/rscale_);
    dy[ii] = dy[ii]*gc;
  }
}

// -----------------------------------------------------------------------------

void ObsLocGC99::print(std::ostream & os) const {
  os << "ObsLocGC99::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
