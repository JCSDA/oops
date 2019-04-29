/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "lorenz95/ObsPreQC.h"

#include <string>

#include "lorenz95/L95Traits.h"
#include "lorenz95/ObsData1D.h"
#include "lorenz95/ObsTableView.h"
#include "oops/interface/ObsFilter.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------
static oops::FilterMaker<L95Traits, oops::ObsFilter<L95Traits, ObsPreQC> > makerPreChk_("PreQC");
// -----------------------------------------------------------------------------

ObsPreQC::ObsPreQC(ObsTableView & obsdb, const eckit::Configuration &,
                   boost::shared_ptr<ObsData1D<int> >,
                   boost::shared_ptr<ObsData1D<float> >) : novars_() {
  if (!obsdb.has("PreQC")) {  // true in MakeObs
    oops::Log::info() << "ObsPreQC::ObsPreQC generating PreQC" << std::endl;
    ObsData1D<int> qc(obsdb, novars_);
    qc.zero();
    qc.save("PreQC");
  }
}

// -----------------------------------------------------------------------------
}  // namespace lorenz95

