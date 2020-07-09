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
    qg_gom_analytic_init_f90(geovals.toFortran(), locs.toFortran(), config_);
}
// -----------------------------------------------------------------------------
}  // namespace qg
