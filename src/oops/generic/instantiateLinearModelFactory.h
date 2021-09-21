/*
 * (C) Copyright 2018-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_INSTANTIATELINEARMODELFACTORY_H_
#define OOPS_GENERIC_INSTANTIATELINEARMODELFACTORY_H_

#include "oops/generic/IdentityLinearModel.h"
#include "oops/generic/LinearModelBase.h"

namespace oops {

template <typename MODEL> void instantiateLinearModelFactory() {
  static LinearModelMaker<MODEL, IdentityLinearModel<MODEL> > makerIdentityLinearModel_("Identity");
}

}  // namespace oops

#endif  // OOPS_GENERIC_INSTANTIATELINEARMODELFACTORY_H_
