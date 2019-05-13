/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/ChangeVarQG.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "model/GeometryQG.h"
#include "model/StateQG.h"
#include "oops/util/Logger.h"

namespace qg {
// -----------------------------------------------------------------------------
ChangeVarQG::ChangeVarQG(const GeometryQG & resol, const eckit::Configuration & conf) {
  oops::Log::trace() << "ChangeVarQG::ChangeVarQG start" << std::endl;
  const VariablesQG vars_in(eckit::LocalConfiguration(conf, "inputVariables"));
  const VariablesQG vars_out(eckit::LocalConfiguration(conf, "outputVariables"));
  qg_change_var_setup_f90(keyConfig_, vars_in.toFortran(), vars_out.toFortran());
  oops::Log::trace() << "ChangeVarQG::ChangeVarQG done" << std::endl;
}
// -----------------------------------------------------------------------------
ChangeVarQG::~ChangeVarQG() {}
// -----------------------------------------------------------------------------
void ChangeVarQG::changeVar(const StateQG & xa, StateQG & xm) const {
  qg_change_var_f90(keyConfig_, xa.fields().toFortran(), xm.fields().toFortran());
}
// -----------------------------------------------------------------------------
void ChangeVarQG::changeVarInverse(const StateQG & xm, StateQG & xa) const {
  qg_change_var_inv_f90(keyConfig_, xm.fields().toFortran(), xa.fields().toFortran());
}
// -----------------------------------------------------------------------------
void ChangeVarQG::print(std::ostream & os) const {
  os << "QG change of variable";
}
// -----------------------------------------------------------------------------
}  // namespace qg


