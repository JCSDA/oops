/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/ChangeVarTLAD.h"

#include <ostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "model/GeometryQG.h"
#include "model/IncrementQG.h"
#include "model/StateQG.h"
#include "oops/util/Logger.h"

namespace qg {
// -----------------------------------------------------------------------------
ChangeVarTLAD::ChangeVarTLAD(const StateQG &, const StateQG &,
                             const GeometryQG & resol, const eckit::Configuration & conf) {
  const eckit::Configuration * configc = &conf;
  qg_setup_f90(&configc, resol.toFortran(), keyConfig_);
  oops::Log::trace() << "ChangeVarTLAD::ChangeVarTLAD constructed" << std::endl;
}
// -----------------------------------------------------------------------------
ChangeVarTLAD::~ChangeVarTLAD() {}
// -----------------------------------------------------------------------------
void ChangeVarTLAD::multiply(const IncrementQG & dxa, IncrementQG & dxm) const {
  dxm = dxa;
  ASSERT(dxm.fields().isForModel(false));
  qg_prepare_integration_tl_f90(keyConfig_, dxm.fields().toFortran());
}
// -----------------------------------------------------------------------------
void ChangeVarTLAD::multiplyInverse(const IncrementQG & dxm, IncrementQG & dxa) const {
  dxa = dxm;
}
// -----------------------------------------------------------------------------
void ChangeVarTLAD::multiplyAD(const IncrementQG & dxm, IncrementQG & dxa) const {
  ASSERT(dxm.fields().isForModel(false));
  qg_prepare_integration_ad_f90(keyConfig_, dxm.fields().toFortran());
  dxa = dxm;
}
// -----------------------------------------------------------------------------
void ChangeVarTLAD::multiplyInverseAD(const IncrementQG & dxa, IncrementQG & dxm) const {
  dxm = dxa;
}
// -----------------------------------------------------------------------------
void ChangeVarTLAD::print(std::ostream & os) const {
  os << "QG linear change of variable";
}
// -----------------------------------------------------------------------------
}  // namespace qg

