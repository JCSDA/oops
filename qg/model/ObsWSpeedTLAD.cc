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

#include <vector>

#include "eckit/config/Configuration.h"
#include "model/GomQG.h"
#include "model/ObsBias.h"
#include "model/ObsBiasIncrement.h"
#include "model/ObsSpaceQG.h"
#include "model/ObsVecQG.h"
#include "model/QgFortran.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
static ObsOpTLADMaker<ObsWSpeedTLAD> makerWSpeedTL_("WSpeed");
// -----------------------------------------------------------------------------

ObsWSpeedTLAD::ObsWSpeedTLAD(const ObsSpaceQG & odb, const eckit::Configuration & config)
  : obsdb_(odb), traj_(), varin_(std::vector<std::string>{"u", "v", "z"})
{
  qg_wspeed_alloctraj_f90(odb.nobs(), traj_);
  oops::Log::trace() << "ObsWSpeedTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsWSpeedTLAD::setTrajectory(const GomQG & gom, const ObsBias &, const QCFlags_ &) {
  qg_wspeed_settraj_f90(obsdb_.toFortran(), gom.toFortran(), traj_);
  oops::Log::trace() << "ObsWSpeedTLAD trajectory was set " << traj_ << std::endl;
}

// -----------------------------------------------------------------------------

void ObsWSpeedTLAD::simulateObsTL(const GomQG & gom, ObsVecQG & ovec,
                                  const ObsBiasIncrement & bias,
                                  const QCFlags_ & qc_flags) const {
  qg_wspeed_equiv_tl_f90(obsdb_.toFortran(), gom.toFortran(), ovec.toFortran(),
                         traj_, bias.wspd());
}

// -----------------------------------------------------------------------------

void ObsWSpeedTLAD::simulateObsAD(GomQG & gom, const ObsVecQG & ovec,
                                  ObsBiasIncrement & bias,
                                  const QCFlags_ & qc_flags) const {
  qg_wspeed_equiv_ad_f90(obsdb_.toFortran(), gom.toFortran(), ovec.toFortran(),
                         traj_, bias.wspd());
}

// -----------------------------------------------------------------------------

void ObsWSpeedTLAD::print(std::ostream & os) const {
  os << "QG wind speed observation operator TL/AD";
}

// -----------------------------------------------------------------------------

}  // namespace qg
