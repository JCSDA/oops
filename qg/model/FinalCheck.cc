/*
 * (C) Crown Copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "model/FinalCheck.h"

#include <string>

#include "model/ObsDataQG.h"
#include "model/ObsSpaceQG.h"
#include "model/QgTraits.h"
#include "oops/interface/ObsFilter.h"

namespace qg {
// -----------------------------------------------------------------------------
static oops::FilterMaker<QgObsTraits,
       oops::ObsFilter<QgObsTraits, FinalCheck> > makerPreChk_("Final Check");
// -----------------------------------------------------------------------------
}  // namespace qg
