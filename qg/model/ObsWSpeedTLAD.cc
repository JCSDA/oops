/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ObsWSpeedTLAD.h"

#include "util/Logger.h"
#include "model/GomQG.h"
#include "model/ObsBias.h"
#include "model/ObsBiasIncrement.h"
#include "model/ObsSpaceQG.h"
#include "model/ObsVecQG.h"
#include "model/QgFortran.h"
#include "model/VariablesQG.h"

using oops::Log;


// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------

ObsWSpeedTLAD::ObsWSpeedTLAD(const ObsSpaceQG & odb, const int & keyOperWspeed)
  : keyOperWspeed_(keyOperWspeed), traj_(), varin_()
{
  int keyVarin;
  qg_obsoper_inputs_f90(keyOperWspeed_, keyVarin);
  varin_.reset(new VariablesQG(keyVarin));
  qg_wspeed_gettraj_f90(keyOperWspeed_, odb.nobs(), traj_.toFortran());
  Log::trace() << "ObsWSpeedTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsWSpeedTLAD::~ObsWSpeedTLAD() {
  Log::trace() << "ObsWSpeedTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsWSpeedTLAD::setTrajectory(const GomQG & gom, const ObsBias &) {
  qg_wspeed_settraj_f90(gom.toFortran(), traj_.toFortran());
  Log::trace() << "ObsWSpeedTLAD trajectory was set " << traj_ << std::endl;
}

// -----------------------------------------------------------------------------

void ObsWSpeedTLAD::obsEquivTL(const GomQG & gom, ObsVecQG & ovec,
                               const ObsBiasIncrement & bias) const {
  Log::debug() << "ObsWSpeedTLAD::obsEquivTL gom " << gom << std::endl;
  Log::debug() << "ObsWSpeedTLAD::obsEquivTL traj " << traj_ << std::endl;
  qg_wspeed_equiv_tl_f90(gom.toFortran(), ovec.toFortran(),
                         traj_.toFortran(), bias.wspd());
  Log::debug() << "ObsWSpeedTLAD::obsEquivTL obsvec " << ovec << std::endl;
}

// -----------------------------------------------------------------------------

void ObsWSpeedTLAD::obsEquivAD(GomQG & gom, const ObsVecQG & ovec,
                               ObsBiasIncrement & bias) const {
  Log::debug() << "ObsWSpeedTLAD::obsEquivAD obsvec " << ovec << std::endl;
  Log::debug() << "ObsWSpeedTLAD::obsEquivTL traj " << traj_ << std::endl;
  qg_wspeed_equiv_ad_f90(gom.toFortran(), ovec.toFortran(),
                         traj_.toFortran(), bias.wspd());
  Log::debug() << "ObsWSpeedTLAD::obsEquivAD gom " << gom << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace qg
