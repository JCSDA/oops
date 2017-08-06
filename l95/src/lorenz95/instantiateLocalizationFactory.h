/*
 * (C) Copyright 2009-2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef LORENZ95_INSTANTIATELOCALIZATIONFACTORY_H_
#define LORENZ95_INSTANTIATELOCALIZATIONFACTORY_H_

#include "oops/interface/LocalizationBase.h"
#include "lorenz95/LocalizationMatrixL95.h"
#include "lorenz95/L95Traits.h"

namespace lorenz95 {

void instantiateLocalizationFactory() {
  static oops::LocalizationMaker<L95Traits, LocalizationMatrixL95> makerWSpeed_("L95");
}

}  // namespace lorenz95

#endif  // LORENZ95_INSTANTIATELOCALIZATIONFACTORY_H_
