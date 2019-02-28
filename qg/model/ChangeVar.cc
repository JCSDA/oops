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
#include "model/StateQG.h"
#include "oops/util/Logger.h"

namespace qg {
// -----------------------------------------------------------------------------
ChangeVar::ChangeVar(const GeometryQG & resol, const eckit::Configuration & conf) {
  oops::Log::trace() << "ChangeVar::ChangeVar start" << std::endl;
  const eckit::Configuration * configc = &conf;
  qg_setup_f90(&configc, resol.toFortran(), keyConfig_);
  oops::Log::trace() << "ChangeVar::ChangeVar done" << std::endl;
}
// -----------------------------------------------------------------------------
ChangeVar::~ChangeVar() {}
// -----------------------------------------------------------------------------
void ChangeVar::changeVar(const StateQG & xa, StateQG & xm) const {
  xm = xa;
  ASSERT(xm.fields().isForModel(false));
  qg_prepare_integration_f90(keyConfig_, xm.fields().toFortran());
  oops::Log::debug() << "ChangeVar::multiply" << xm << std::endl;
}
// -----------------------------------------------------------------------------
void ChangeVar::changeVarInverse(const StateQG & xm, StateQG & xa) const {
  xa = xm;
}
// -----------------------------------------------------------------------------
void ChangeVar::print(std::ostream & os) const {
  os << "QG change of variable";
}
// -----------------------------------------------------------------------------
}  // namespace qg

