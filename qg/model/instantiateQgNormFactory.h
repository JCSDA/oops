/*
 * (C) Crown Copyright 2023, the Met Office
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef QG_MODEL_INSTANTIATEQGNORMFACTORY_H_
#define QG_MODEL_INSTANTIATEQGNORMFACTORY_H_


#include "oops/base/NormBase.h"
#include "oops/generic/instantiateNormFactory.h"
#include "oops/generic/NormScalar.h"
#include "oops/qg/QgTraits.h"

namespace qg {

void instantiateQgNormFactory() {
  oops::instantiateNormFactory<qg::QgTraits>();
}

}  // namespace qg

#endif  // QG_MODEL_INSTANTIATEQGNORMFACTORY_H_
