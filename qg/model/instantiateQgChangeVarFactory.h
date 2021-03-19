/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_INSTANTIATEQGCHANGEVARFACTORY_H_
#define QG_MODEL_INSTANTIATEQGCHANGEVARFACTORY_H_

#include "oops/interface/LinearVariableChange.h"
#include "oops/interface/VariableChange.h"

#include "oops/qg/ChangeVarQG.h"
#include "oops/qg/ChangeVarTLADQG.h"
#include "oops/qg/QgTraits.h"

namespace qg {

void instantiateQgChangeVarFactory() {
  static oops::VariableChangeMaker<QgTraits, ChangeVarQG> makerchangevar_("ChVarQG");
  static oops::VariableChangeMaker<QgTraits, ChangeVarQG> makerdefchavar_("default");

  static oops::LinearVariableChangeMaker<qg::QgTraits,
                                   oops::LinearVariableChange<qg::QgTraits, qg::ChangeVarTLADQG> >
               makerChLinVarQG_("ChVarQG");
}

}  // namespace qg

#endif  // QG_MODEL_INSTANTIATEQGCHANGEVARFACTORY_H_
