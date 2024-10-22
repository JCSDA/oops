/*
 * (C) Copyright 2024- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/ObsFilter.h"

#include "eckit/config/Configuration.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
ObsFilter::ObsFilter(const ObsSpaceQG &, const eckit::Configuration & conf,
             std::shared_ptr<ObsDataQG<int> >, std::shared_ptr<ObsDataQG<float> >,
             const int iteration)
  : novars_(), config_(conf.getSubConfiguration("obs filtering")),
    saveGeoVaLs_(config_.getBool("save geovals", false))
{}

// -----------------------------------------------------------------------------
void ObsFilter::priorFilter(const GomQG & gv) {
  if (saveGeoVaLs_) {
    gv.write(config_);
  }
}

}  // namespace qg

