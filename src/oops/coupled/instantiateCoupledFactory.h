/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "oops/generic/ModelBase.h"

#include "oops/coupled/ModelCoupled.h"
#include "oops/coupled/TraitCoupled.h"

namespace oops {

template <typename MODEL1, typename MODEL2>
void instantiateCoupledFactory() {
  static interface::ModelMaker<TraitCoupled<MODEL1, MODEL2>,
                               ModelCoupled<MODEL1, MODEL2>
                              > makerModelCoupled_("Coupled");
}

}  // namespace oops

