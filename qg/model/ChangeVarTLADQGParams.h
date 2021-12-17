/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QG_MODEL_CHANGEVARTLADQGPARAMS_H_
#define QG_MODEL_CHANGEVARTLADQGPARAMS_H_

#include <ostream>
#include <string>

#include "oops/base/LinearVariableChangeParametersBase.h"

namespace qg {

// -------------------------------------------------------------------------------------------------
/// No change of linear variable parameters

class ChangeVarTLADQGParams : public oops::LinearVariableChangeParametersBase {
  OOPS_CONCRETE_PARAMETERS(ChangeVarTLADQGParams, LinearVariableChangeParametersBase)
 public:
  // Linear variable change parameters would go here.
};

// -------------------------------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_CHANGEVARTLADQGPARAMS_H_
