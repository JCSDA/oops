/*
 * (C) Copyright 2017-2021 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef LORENZ95_INSTANTIATELOCALIZATIONFACTORY_H_
#define LORENZ95_INSTANTIATELOCALIZATIONFACTORY_H_

#include "oops/generic/LocalizationBase.h"
#include "oops/interface/Localization.h"

namespace lorenz95 {

template<typename L95MODEL>
void instantiateLocalizationFactory() {
  static oops::LocalizationMaker<L95MODEL, oops::interface::Localization<L95MODEL>>
      makerL95loc_("L95");
}

}  // namespace lorenz95

#endif  // LORENZ95_INSTANTIATELOCALIZATIONFACTORY_H_
