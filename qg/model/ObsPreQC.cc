/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "model/ObsPreQC.h"

#include <string>

#include "eckit/config/Configuration.h"

#include "model/ObsSpaceQG.h"
#include "model/ObsVecQG.h"
#include "model/QgTraits.h"
#include "oops/interface/ObsFilter.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
static oops::FilterMaker<QgTraits, oops::ObsFilter<QgTraits, ObsPreQC> >
  makerPreChk_("PreQC");
// -----------------------------------------------------------------------------

ObsPreQC::ObsPreQC(ObsSpaceQG & obsdb, const eckit::Configuration & config) {
  const std::string qcname(config.getString("QCname"));
  const oops::Variables var(config.getStringVector("observed"));
  ObsVecQG qc(obsdb, var);
  qc.save(qcname);
}

// -----------------------------------------------------------------------------

}  // namespace qg
