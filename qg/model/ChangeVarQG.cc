/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/ChangeVarQG.h"

#include <ostream>
#include <string>

#include "oops/util/Logger.h"

#include "model/QgFortran.h"
#include "model/StateQG.h"

namespace qg {
// -----------------------------------------------------------------------------
ChangeVarQG::ChangeVarQG(const GeometryQG &, const eckit::Configuration &) {}
// -----------------------------------------------------------------------------
ChangeVarQG::~ChangeVarQG() {}
// -----------------------------------------------------------------------------
void ChangeVarQG::changeVar(const StateQG & xa, StateQG & xm) const {
  qg_change_var_f90(xa.fields().toFortran(), xm.fields().toFortran());
}
// -----------------------------------------------------------------------------
void ChangeVarQG::changeVarInverse(const StateQG & xm, StateQG & xa) const {
  qg_change_var_f90(xm.fields().toFortran(), xa.fields().toFortran());
}
// -----------------------------------------------------------------------------
void ChangeVarQG::print(std::ostream & os) const {
  os << "QG change of variable";
}
// -----------------------------------------------------------------------------
}  // namespace qg


