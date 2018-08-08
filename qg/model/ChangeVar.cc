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
// static oops::VariableChangeMaker<QgTraits, ChangeVar> makerChVarQG("QgCV");
// -----------------------------------------------------------------------------
ChangeVar::ChangeVar(const eckit::Configuration &) {}
// -----------------------------------------------------------------------------
ChangeVar::~ChangeVar() {}
// -----------------------------------------------------------------------------
void ChangeVar::linearize(const StateQG &, const GeometryQG &) {}
// -----------------------------------------------------------------------------
void ChangeVar::transform(const IncrementQG & dxa, IncrementQG & dxm) const {
  dxm = dxa;
}
// -----------------------------------------------------------------------------
void ChangeVar::transformInverse(const IncrementQG & dxm, IncrementQG & dxa) const {
  dxa = dxm;
}
// -----------------------------------------------------------------------------
void ChangeVar::transformAdjoint(const IncrementQG & dxm, IncrementQG & dxa) const {
  dxa = dxm;
}
// -----------------------------------------------------------------------------
void ChangeVar::transformInverseAdjoint(const IncrementQG & dxa, IncrementQG & dxm) const {
  dxm = dxa;
}
// -----------------------------------------------------------------------------
void ChangeVar::print(std::ostream & os) const {
  os << "QG change variable";
}
// -----------------------------------------------------------------------------
}  // namespace qg

