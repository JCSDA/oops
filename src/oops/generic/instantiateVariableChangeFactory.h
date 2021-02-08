/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_GENERIC_INSTANTIATEVARIABLECHANGEFACTORY_H_
#define OOPS_GENERIC_INSTANTIATEVARIABLECHANGEFACTORY_H_

#include "oops/base/LinearVariableChangeBase.h"
#include "oops/base/VariableChangeBase.h"
#include "oops/generic/IdLinearVariableChange.h"
#include "oops/generic/IdVariableChange.h"

namespace oops {

template <typename MODEL>
void instantiateVariableChangeFactory() {
// Nonlinear change of variables
  static GenericVariableChangeMaker<MODEL, IdVariableChange<MODEL> > makerId_("Identity");

// Linear change of variables
  static LinearVariableChangeMaker<MODEL, IdLinearVariableChange<MODEL> > makerIdLin_("Identity");
}

}  // namespace oops

#endif  // OOPS_GENERIC_INSTANTIATEVARIABLECHANGEFACTORY_H_
