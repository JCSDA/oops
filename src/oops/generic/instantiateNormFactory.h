/*
 * (C) Crown Copyright 2023, the Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_INSTANTIATENORMFACTORY_H_
#define OOPS_GENERIC_INSTANTIATENORMFACTORY_H_

#include "oops/base/NormBase.h"
#include "oops/generic/NormScalar.h"

namespace oops {

template <typename MODEL> void instantiateNormFactory() {
  static NormMaker<MODEL, NormScalar<MODEL> >     makerNormScalar_("Scalar");
}

}  // namespace oops

#endif  // OOPS_GENERIC_INSTANTIATENORMFACTORY_H_
