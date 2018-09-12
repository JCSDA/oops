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
#include "oops/util/Logger.h"

namespace qg {
// -----------------------------------------------------------------------------
ChangeVar::ChangeVar(const StateQG &, const StateQG &,
                     const GeometryQG & resol, const eckit::Configuration & conf) {
  oops::Log::trace() << "ChangeVar::ChangeVar start" << std::endl;
  const eckit::Configuration * configc = &conf;
  qg_setup_f90(&configc, resol.toFortran(), keyConfig_);
  oops::Log::trace() << "ChangeVar::ChangeVar done" << std::endl;
}
// -----------------------------------------------------------------------------
ChangeVar::~ChangeVar() {}
// -----------------------------------------------------------------------------
void ChangeVar::multiply(const IncrementQG & dxa, IncrementQG & dxm) const {
  dxm = dxa;
//  ASSERT(dxm.fields().isForModel(false));
//  qg_prepare_integration_tl_f90(keyConfig_, dxm.fields().toFortran());
  oops::Log::debug() << "ChangeVar::multiply" << dxm << std::endl;
}
// -----------------------------------------------------------------------------
void ChangeVar::multiplyInverse(const IncrementQG & dxm, IncrementQG & dxa) const {
  dxa = dxm;
}
// -----------------------------------------------------------------------------
void ChangeVar::multiplyAD(const IncrementQG & dxm, IncrementQG & dxa) const {
//  ASSERT(dxm.fields().isForModel(false));
//  qg_prepare_integration_ad_f90(keyConfig_, dxm.fields().toFortran());
  dxa = dxm;
  oops::Log::debug() << "ChangeVar::multiplyAD" << dxa << std::endl;
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

