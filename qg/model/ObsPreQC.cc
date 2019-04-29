/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "model/ObsPreQC.h"

#include <string>

#include "model/ObsDataQG.h"
#include "model/ObsSpaceQG.h"
#include "model/QgTraits.h"
#include "oops/interface/ObsFilter.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
static oops::FilterMaker<QgTraits, oops::ObsFilter<QgTraits, ObsPreQC> > makerPreChk_("PreQC");
// -----------------------------------------------------------------------------

ObsPreQC::ObsPreQC(ObsSpaceQG & obsdb, const eckit::Configuration &,
                   boost::shared_ptr<ObsDataQG<int> >,
                   boost::shared_ptr<ObsDataQG<float> >) : novars_() {
  if (!obsdb.has("PreQC")) {  // true in MakeObs
    oops::Log::info() << "ObsPreQC::ObsPreQC generating PreQC" << std::endl;
    ObsDataQG<int> qc(obsdb, novars_);
    qc.zero();
    qc.save("PreQC");
  }
}

// -----------------------------------------------------------------------------

}  // namespace qg
