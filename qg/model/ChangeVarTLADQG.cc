/*
 * (C) Copyright 2017-2021  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/ChangeVarTLADQG.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "model/GeometryQG.h"
#include "model/IncrementQG.h"
#include "model/QgFortran.h"
#include "model/StateQG.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace qg {
// -----------------------------------------------------------------------------
ChangeVarTLADQG::ChangeVarTLADQG(const GeometryQG &, const Parameters_ &) {}
// -----------------------------------------------------------------------------
ChangeVarTLADQG::~ChangeVarTLADQG() {}
// -----------------------------------------------------------------------------
void ChangeVarTLADQG::changeVarTL(IncrementQG & dx, const oops::Variables & vars) const {
  qg_change_var_tl_f90(dx.fields().toFortran(), vars);
}
// -----------------------------------------------------------------------------
void ChangeVarTLADQG::changeVarInverseTL(IncrementQG & dx, const oops::Variables & vars) const {
  qg_change_var_tl_f90(dx.fields().toFortran(), vars);
}
// -----------------------------------------------------------------------------
void ChangeVarTLADQG::changeVarAD(IncrementQG & dx, const oops::Variables & vars) const {
  qg_change_var_ad_f90(dx.fields().toFortran(), vars);
}
// -----------------------------------------------------------------------------
void ChangeVarTLADQG::changeVarInverseAD(IncrementQG & dx, const oops::Variables & vars) const {
  qg_change_var_ad_f90(dx.fields().toFortran(), vars);
}
// -----------------------------------------------------------------------------
void ChangeVarTLADQG::changeVarTraj(const StateQG & firstGuess, const oops::Variables & vars) {
  // No QG trajectory used. No fortran to call here.
}
// -----------------------------------------------------------------------------
void ChangeVarTLADQG::print(std::ostream & os) const {
  os << "QG linear change of variable";
}
// -----------------------------------------------------------------------------
}  // namespace qg

