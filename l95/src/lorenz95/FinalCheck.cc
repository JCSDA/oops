/*
 * (C) Crown Copyright 2021, Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "lorenz95/FinalCheck.h"

#include <string>

#include "lorenz95/L95Traits.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------
static oops::interface::FilterMaker<L95ObsTraits, FinalCheck> makerQCm_("Final Check");
// -----------------------------------------------------------------------------
}  // namespace lorenz95

