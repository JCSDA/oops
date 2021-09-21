/*
 * (C) Copyright 2018-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_INSTANTIATEMODELFACTORY_H_
#define OOPS_GENERIC_INSTANTIATEMODELFACTORY_H_

#include "oops/generic/IdentityModel.h"
#include "oops/generic/ModelBase.h"
#include "oops/generic/PseudoModel.h"

namespace oops {

template <typename MODEL> void instantiateModelFactory() {
  static ModelMaker<MODEL, IdentityModel<MODEL> > makerIdentityModel_("Identity");
  static ModelMaker<MODEL, PseudoModel<MODEL> > makerPseudoModel_("PseudoModel");
}

}  // namespace oops

#endif  // OOPS_GENERIC_INSTANTIATEMODELFACTORY_H_
