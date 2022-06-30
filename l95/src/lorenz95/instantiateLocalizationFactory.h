/*
 * (C) Copyright 2017-2021 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef LORENZ95_INSTANTIATELOCALIZATIONFACTORY_H_
#define LORENZ95_INSTANTIATELOCALIZATIONFACTORY_H_

#include "lorenz95/L95Traits.h"
#include "lorenz95/LocalizationMatrixL95.h"
#include "oops/interface/LocalizationBase.h"

namespace lorenz95 {

void instantiateLocalizationFactory() {
  static oops::interface::LocalizationMaker<L95Traits, LocalizationMatrixL95> makerL95_("L95");
}

}  // namespace lorenz95

#endif  // LORENZ95_INSTANTIATELOCALIZATIONFACTORY_H_
