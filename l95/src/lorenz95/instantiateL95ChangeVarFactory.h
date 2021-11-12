/*
 * (C) Copyright 2017-2020 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef LORENZ95_INSTANTIATEL95CHANGEVARFACTORY_H_
#define LORENZ95_INSTANTIATEL95CHANGEVARFACTORY_H_

#include "lorenz95/L95Traits.h"
#include "oops/base/VariableChangeBase.h"
#include "oops/generic/IdLinearVariableChange.h"
#include "oops/generic/IdVariableChange.h"
#include "oops/interface/LinearVariableChange.h"

namespace lorenz95 {

void instantiateL95ChangeVarFactory() {
  static oops::GenericVariableChangeMaker<L95Traits, oops::IdVariableChange<L95Traits> >
           makerL95_("default");

  static oops::LinearVariableChangeMaker<L95Traits, oops::IdLinearVariableChange<L95Traits> >
           makerL95LinVarDef_("default");
}

}  // namespace lorenz95

#endif  // LORENZ95_INSTANTIATEL95CHANGEVARFACTORY_H_
