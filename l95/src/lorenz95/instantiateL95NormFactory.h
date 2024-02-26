/*
 * (C) Crown Copyright 2023, the Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef LORENZ95_INSTANTIATEL95NORMFACTORY_H_
#define LORENZ95_INSTANTIATEL95NORMFACTORY_H_

#include "lorenz95/L95Traits.h"
#include "oops/base/Norm.h"
#include "oops/generic/instantiateNormFactory.h"
#include "oops/generic/NormScalar.h"

namespace lorenz95 {

void instantiateL95NormFactory() {
  oops::instantiateNormFactory<L95Traits>();
}

}  // namespace lorenz95

#endif  // LORENZ95_INSTANTIATEL95NORMFACTORY_H_
