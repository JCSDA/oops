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

#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "model/GeometryQG.h"
#include "model/IncrementQG.h"
#include "model/ModelBiasIncrement.h"
#include "model/ModelQG.h"
#include "model/QgFortran.h"
#include "model/QgTraits.h"
#include "model/StateQG.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace qg {
// -----------------------------------------------------------------------------
static oops::LinearModelMaker<QgTraits, TlmQG> makerQGTLM_("QgTLM");
// -----------------------------------------------------------------------------
TlmQG::TlmQG(const GeometryQG & resol, const eckit::Configuration & tlConf)
  : keyConfig_(0), tstep_(), resol_(resol), traj_(),
    lrmodel_(resol_, eckit::LocalConfiguration(tlConf, "trajectory")),
    linvars_(std::vector<std::string>{"x", "q", "u", "v"})
{
  tstep_ = util::Duration(tlConf.getString("tstep"));

  const eckit::Configuration * configc = &tlConf;
  qg_setup_f90(&configc, resol_.toFortran(), keyConfig_);

  oops::Log::trace() << "TlmQG created" << std::endl;
}
// -----------------------------------------------------------------------------
TlmQG::~TlmQG() {
  qg_delete_f90(keyConfig_);
  for (trajIter jtra = traj_.begin(); jtra != traj_.end(); ++jtra) {
    qg_wipe_traj_f90(jtra->second);
  }
  oops::Log::trace() << "TlmQG destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void TlmQG::setTrajectory(const StateQG & xx, StateQG & xlr, const ModelBias & bias) {
// StateQG xlr(resol_, xx);
  xlr.changeResolution(xx);
  int ftraj = lrmodel_.saveTrajectory(xlr, bias);
  traj_[xx.validTime()] = ftraj;

// should be in print method
  std::vector<double> zstat(15);
  qg_traj_minmaxrms_f90(ftraj, zstat[0]);
  oops::Log::debug() << "TlmQG trajectory at time " << xx.validTime() << std::endl;
  for (unsigned int jj = 0; jj < 5; ++jj) {
    oops::Log::debug() << "  Min=" << zstat[3*jj] << ", Max=" << zstat[3*jj+1]
                       << ", RMS=" << zstat[3*jj+2] << std::endl;
  }
// should be in print method
}
// -----------------------------------------------------------------------------
void TlmQG::initializeTL(IncrementQG & dx) const {
  ASSERT(dx.fields().isForModel(false));
  qg_prepare_integration_tl_f90(keyConfig_, dx.fields().toFortran());
  oops::Log::debug() << "TlmQG::initializeTL" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmQG::stepTL(IncrementQG & dx, const ModelBiasIncrement &) const {
  trajICst itra = traj_.find(dx.validTime());
  if (itra == traj_.end()) {
    oops::Log::error() << "TlmQG: trajectory not available at time " << dx.validTime() << std::endl;
    ABORT("TlmQG: trajectory not available");
  }
  ASSERT(dx.fields().isForModel(false));
  oops::Log::debug() << "TlmQG::stepTL fields in" << dx.fields() << std::endl;
  qg_propagate_tl_f90(keyConfig_, dx.fields().toFortran(), itra->second);
  oops::Log::debug() << "TlmQG::stepTL fields out" << dx.fields() << std::endl;
  dx.validTime() += tstep_;
}
// -----------------------------------------------------------------------------
void TlmQG::finalizeTL(IncrementQG & dx) const {
  oops::Log::debug() << "TlmQG::finalizeTL" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmQG::initializeAD(IncrementQG & dx) const {
  ASSERT(dx.fields().isForModel(false));
  oops::Log::debug() << "TlmQG::initializeAD" << dx.fields() << std::endl;
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
  oops::Log::debug() << "TlmQG::stepAD fields in" << dx.fields() << std::endl;
  qg_propagate_ad_f90(keyConfig_, dx.fields().toFortran(), itra->second);
  oops::Log::debug() << "TlmQG::stepAD fields out" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmQG::finalizeAD(IncrementQG & dx) const {
  qg_prepare_integration_ad_f90(keyConfig_, dx.fields().toFortran());
  oops::Log::debug() << "TlmQG::finalizeAD" << dx.fields() << std::endl;
}
// -----------------------------------------------------------------------------
void TlmQG::print(std::ostream & os) const {
  os << "QG TLM Trajectory, nstep=" << traj_.size() << std::endl;
  typedef std::map< util::DateTime, int >::const_iterator trajICst;
  if (traj_.size() > 0) {
    os << "QG TLM Trajectory: times are:";
    for (trajICst jtra = traj_.begin(); jtra != traj_.end(); ++jtra) {
      os << "  " << jtra->first;
    }
  }
}
// -----------------------------------------------------------------------------
}  // namespace qg
