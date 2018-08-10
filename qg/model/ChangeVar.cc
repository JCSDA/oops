/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/ChangeVar.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "model/GeometryQG.h"
#include "model/IncrementQG.h"
#include "model/StateQG.h"

namespace qg {
// -----------------------------------------------------------------------------
ChangeVar::ChangeVar(const eckit::Configuration &) {}
// -----------------------------------------------------------------------------
ChangeVar::~ChangeVar() {}
// -----------------------------------------------------------------------------
void ChangeVar::linearize(const StateQG &, const GeometryQG &) {}
// -----------------------------------------------------------------------------
void ChangeVar::multiply(const IncrementQG & dxa, IncrementQG & dxm) const {
  dxm = dxa;
}
// -----------------------------------------------------------------------------
void ChangeVar::multiplyInverse(const IncrementQG & dxm, IncrementQG & dxa) const {
  dxa = dxm;
}
// -----------------------------------------------------------------------------
void ChangeVar::multiplyAD(const IncrementQG & dxm, IncrementQG & dxa) const {
  dxa = dxm;
}
// -----------------------------------------------------------------------------
void ChangeVar::multiplyInverseAD(const IncrementQG & dxa, IncrementQG & dxm) const {
  dxm = dxa;
}
// -----------------------------------------------------------------------------
void ChangeVar::print(std::ostream & os) const {
  os << "QG change variable";
}
// -----------------------------------------------------------------------------
}  // namespace qg

