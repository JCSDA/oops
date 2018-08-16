/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <vector>

#include "model/TlmIdQG.h"

#include "eckit/config/LocalConfiguration.h"
#include "model/GeometryQG.h"
#include "model/IncrementQG.h"
#include "model/ModelBiasIncrement.h"
#include "model/QgFortran.h"
#include "model/QgTraits.h"
#include "model/StateQG.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace qg {
// -----------------------------------------------------------------------------
static oops::LinearModelMaker<QgTraits, TlmIdQG> makerQGIdTLM_("QgIdTLM");
// -----------------------------------------------------------------------------
TlmIdQG::TlmIdQG(const GeometryQG & resol, const eckit::Configuration & tlConf)
  : keyConfig_(0), tstep_(), resol_(resol),
    linvars_(std::vector<std::string>{"x", "q", "u", "v"})
{
  tstep_ = util::Duration(tlConf.getString("tstep"));

  const eckit::Configuration * configc = &tlConf;
  qg_setup_f90(&configc, resol_.toFortran(), keyConfig_);

  oops::Log::trace() << "TlmIdQG created" << std::endl;
}
// -----------------------------------------------------------------------------
TlmIdQG::~TlmIdQG() {
  qg_delete_f90(keyConfig_);
  oops::Log::trace() << "TlmIdQG destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdQG::setTrajectory(const StateQG &, StateQG &, const ModelBias &) {}
// -----------------------------------------------------------------------------
void TlmIdQG::initializeTL(IncrementQG & dx) const {
  dx.activateModel();
  ASSERT(dx.fields().isForModel(false));
  qg_prepare_integration_tl_f90(keyConfig_, dx.fields().toFortran());
  oops::Log::debug() << "TlmIdQG::initializeTL" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdQG::stepTL(IncrementQG & dx, const ModelBiasIncrement &) const {
  dx.updateTime(tstep_);
}
// -----------------------------------------------------------------------------
void TlmIdQG::finalizeTL(IncrementQG & dx) const {
  dx.deactivateModel();
  oops::Log::debug() << "TlmIdQG::finalizeTL" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdQG::initializeAD(IncrementQG & dx) const {
  dx.activateModel();
  ASSERT(dx.fields().isForModel(false));
  oops::Log::debug() << "TlmIdQG::initializeAD" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdQG::stepAD(IncrementQG & dx, ModelBiasIncrement &) const {
  dx.updateTime(-tstep_);
}
// -----------------------------------------------------------------------------
void TlmIdQG::finalizeAD(IncrementQG & dx) const {
  qg_prepare_integration_ad_f90(keyConfig_, dx.fields().toFortran());
  dx.deactivateModel();
  oops::Log::debug() << "TlmIdQG::finalizeAD" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmIdQG::print(std::ostream & os) const {
  os << "QG IdTLM" << std::endl;
}
// -----------------------------------------------------------------------------
}  // namespace qg
