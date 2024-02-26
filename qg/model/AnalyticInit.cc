/*
 * (C) Copyright 2020-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/AnalyticInit.h"

#include "model/GomQG.h"
#include "model/LocationsQG.h"

namespace qg {

// Ideally, we would have two separate analytic init classes; one for baroclinic
// instability, one for large vortices, but for now we'll use the same class to
// do both.
static oops::AnalyticInitMaker<QgObsTraits, AnalyticInit> makerAnalytic1_("baroclinic-instability");
static oops::AnalyticInitMaker<QgObsTraits, AnalyticInit> makerAnalytic2_("large-vortices");

// -----------------------------------------------------------------------------
AnalyticInit::AnalyticInit(const eckit::Configuration & config): config_(config)
{ }
// -----------------------------------------------------------------------------
/*! \brief GeoVaLs Analytic Initialization
 *
 * \details This qg::AnalyticInit constructor was introduced in May, 2018 for use with
 * the interpolation test.
 *
 */
void AnalyticInit::fillGeoVaLs(const LocationsQG & locs,
                               GomQG & geovals) const {
  // Optionally replace values with analytic init
  if (config_.has("analytic_init"))
    qg_gom_analytic_init_f90(geovals.toFortran(), locs, config_);
}
// -----------------------------------------------------------------------------
}  // namespace qg
