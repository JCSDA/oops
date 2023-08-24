/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/TlmQG.h"

#include <iomanip>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

#include "model/GeometryQG.h"
#include "model/IncrementQG.h"
#include "model/ModelBias.h"
#include "model/ModelBiasIncrement.h"
#include "model/ModelQG.h"
#include "model/QgFortran.h"
#include "model/QgTraits.h"
#include "model/StateQG.h"

namespace qg {
// -----------------------------------------------------------------------------
static oops::interface::LinearModelMaker<QgTraits, TlmQG> makerQGTLM_("QgTLM");
// -----------------------------------------------------------------------------
TlmQG::TlmQG(const GeometryQG & resol, const eckit::Configuration & tlConf)
  : keyConfig_(0), tstep_(util::Duration(tlConf.getString("tstep"))), resol_(resol), traj_(),
    lrmodel_(resol_, eckit::LocalConfiguration(tlConf, "trajectory")),
    linvars_({"x"})
{
  if (tlConf.has("tlm variables")) linvars_ = oops::Variables(tlConf, "tlm variables");
  qg_model_setup_f90(keyConfig_, tlConf);

  oops::Log::trace() << "TlmQG created" << std::endl;
}
// -----------------------------------------------------------------------------
TlmQG::~TlmQG() {
  qg_model_delete_f90(keyConfig_);
  for (trajIter jtra = traj_.begin(); jtra != traj_.end(); ++jtra) {
    qg_fields_delete_f90(jtra->second);
  }
  oops::Log::trace() << "TlmQG destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmQG::setTrajectory(const StateQG & xx, StateQG & xlr, const ModelBias & bias) {
// StateQG xlr(resol_, xx);
  xlr.changeResolution(xx);
  int ftraj = lrmodel_.saveTrajectory(xlr, bias);
  traj_[xx.validTime()] = ftraj;
}
// -----------------------------------------------------------------------------
void TlmQG::initializeTL(IncrementQG & dx) const {
  ASSERT(dx.fields().isForModel(false));
}
// -----------------------------------------------------------------------------
void TlmQG::stepTL(IncrementQG & dx, const ModelBiasIncrement &) const {
  trajICst itra = traj_.find(dx.validTime());
  if (itra == traj_.end()) {
    oops::Log::error() << "TlmQG: trajectory not available at time " << dx.validTime() << std::endl;
    ABORT("TlmQG: trajectory not available");
  }
  ASSERT(dx.fields().isForModel(false));
  qg_model_propagate_tl_f90(keyConfig_, itra->second, dx.fields().toFortran());
  dx.validTime() += tstep_;
}
// -----------------------------------------------------------------------------
void TlmQG::finalizeTL(IncrementQG & dx) const {}
// -----------------------------------------------------------------------------
void TlmQG::initializeAD(IncrementQG & dx) const {
  ASSERT(dx.fields().isForModel(false));
}
// -----------------------------------------------------------------------------
void TlmQG::stepAD(IncrementQG & dx, ModelBiasIncrement &) const {
  dx.validTime() -= tstep_;
  trajICst itra = traj_.find(dx.validTime());
  if (itra == traj_.end()) {
    oops::Log::error() << "TlmQG: trajectory not available at time " << dx.validTime() << std::endl;
    ABORT("TlmQG: trajectory not available");
  }
  ASSERT(dx.fields().isForModel(false));
  qg_model_propagate_ad_f90(keyConfig_, itra->second, dx.fields().toFortran());
}
// -----------------------------------------------------------------------------
void TlmQG::finalizeAD(IncrementQG & dx) const {}
// -----------------------------------------------------------------------------
void TlmQG::print(std::ostream & os) const {
  os << "QG TLM Trajectory, nstep=" << traj_.size() << std::endl;
  typedef std::map< util::DateTime, int >::const_iterator trajICst;
  if (traj_.size() > 0) {
    os << "QG TLM Trajectory: times are:";
    for (trajICst jtra = traj_.begin(); jtra != traj_.end(); ++jtra) {
      // Time
      os << "  " << jtra->first;
      os << "  " << jtra->second;
    }
  }
}
// -----------------------------------------------------------------------------
}  // namespace qg
