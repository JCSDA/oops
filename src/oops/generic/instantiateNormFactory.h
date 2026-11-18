/*
 * (C) Crown Copyright 2026, Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#pragma once

#include "oops/base/NormBase.h"
#include "oops/generic/NormScalar.h"

namespace oops {

template <typename MODEL> void instantiateNormFactory() {
  static typename NormFactory<MODEL>::template Maker<NormScalar<MODEL>> makerNormScalar_("Scalar");

  MODEL::registerNorms();
}

}  // namespace oops
