/*
 * (C) Copyright 2017-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/ChangeVarQG.h"

#include <ostream>
#include <string>

#include "model/QgFortran.h"
#include "model/StateQG.h"

namespace qg {
// -----------------------------------------------------------------------------
ChangeVarQG::ChangeVarQG(const eckit::Configuration &, const GeometryQG &) {}
// -----------------------------------------------------------------------------
ChangeVarQG::~ChangeVarQG() {}
// -----------------------------------------------------------------------------
void ChangeVarQG::changeVar(StateQG & xx, const oops::Variables & vars) const {
  qg_change_var_f90(xx.fields().toFortran(), vars);
  xx.fields().variables() = vars;
}
// -----------------------------------------------------------------------------
void ChangeVarQG::changeVarInverse(StateQG & xx, const oops::Variables & vars) const {
  qg_change_var_f90(xx.fields().toFortran(), vars);
  xx.fields().variables() = vars;
}
// -----------------------------------------------------------------------------
void ChangeVarQG::print(std::ostream & os) const {
  os << "QG change of variable";
}
// -----------------------------------------------------------------------------
}  // namespace qg
