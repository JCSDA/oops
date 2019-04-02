/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "lorenz95/ObsPreQC.h"

#include <string>

#include "eckit/config/Configuration.h"

#include "lorenz95/L95Traits.h"
#include "lorenz95/ObsTableView.h"
#include "lorenz95/ObsVec1D.h"
#include "oops/interface/ObsFilter.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------
static oops::FilterMaker<L95Traits, oops::ObsFilter<L95Traits, ObsPreQC> >
  makerPreChk_("PreQC");
// -----------------------------------------------------------------------------

ObsPreQC::ObsPreQC(ObsTableView & obsdb, const eckit::Configuration & config)
     : novars_() {
  const std::string qcname(config.getString("QCname"));
  const oops::Variables var(config.getStringVector("observed"));
  ObsVec1D qc(obsdb, var);
  qc.save(qcname);
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
