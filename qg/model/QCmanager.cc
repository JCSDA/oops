/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "model/QCmanager.h"

namespace qg {
// -----------------------------------------------------------------------------
static oops::interface::FilterMaker<QgObsTraits, QCmanager> makerPreChk_("QCmanager");
// -----------------------------------------------------------------------------
}  // namespace qg
