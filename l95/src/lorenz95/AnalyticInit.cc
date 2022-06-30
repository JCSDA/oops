/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "lorenz95/AnalyticInit.h"

#include <cmath>

#include "lorenz95/GomL95.h"
#include "lorenz95/LocsL95.h"

#include "oops/util/Logger.h"

namespace lorenz95 {

static oops::AnalyticInitMaker<L95ObsTraits, AnalyticInit> makerAnalytic_("sinus");

// -----------------------------------------------------------------------------
AnalyticInit::AnalyticInit(const Parameters_ & options): options_(options)
{ }
// -----------------------------------------------------------------------------
/*! GomL95 analytic initialization */
void AnalyticInit::fillGeoVaLs(const LocsL95 & locs,
                               GomL95 & geovals) const {
  oops::Log::trace() << "AnalyticInit::fillGeoVaLs " << std::endl;
  const size_t nlocs = geovals.size();
  // analytic init for testing interpolation
  // mean value; default is zero.
  for (size_t jj = 0; jj < nlocs; ++jj) geovals[jj] = options_.mean;
  // add a sinus function if specified in yaml
  if (options_.sinus.value() != boost::none) {
    const double zz = *options_.sinus.value();
    for (size_t jj = 0; jj < nlocs; ++jj)
      geovals[jj] += zz * std::sin(2.0*M_PI*locs[jj]);
  } else if (options_.mean == 0.0) {
    oops::Log::warning() << "Using default analytic init (all zeros)" << std::endl;
  }

  oops::Log::trace() << "GomL95::GomL95 analytic init finished" << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace lorenz95
