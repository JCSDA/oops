/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ObsWindTLAD.h"

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
static ObsOpTLADMaker<ObsWindTLAD> makerWindTL_("Wind");
// -----------------------------------------------------------------------------

ObsWindTLAD::ObsWindTLAD(const ObsSpaceQG & odb, const eckit::Configuration & config)
  : obsdb_(odb), varin_(std::vector<std::string>{"u", "v", "z"})
{
  oops::Log::trace() << "ObsWindTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsWindTLAD::setTrajectory(const GomQG &, const ObsBias &, const QCFlags_ &) {}

// -----------------------------------------------------------------------------

void ObsWindTLAD::simulateObsTL(const GomQG & gom, ObsVecQG & ovec,
                                const ObsBiasIncrement & bias,
                                const QCFlags_ & qc_flags) const {
  qg_wind_equiv_tl_f90(obsdb_.toFortran(), gom.toFortran(), ovec.toFortran(), bias.wind());
}

// -----------------------------------------------------------------------------

void ObsWindTLAD::simulateObsAD(GomQG & gom, const ObsVecQG & ovec,
                                ObsBiasIncrement & bias,
                                const QCFlags_ & qc_flags) const {
  qg_wind_equiv_ad_f90(obsdb_.toFortran(), gom.toFortran(), ovec.toFortran(), bias.wind());
}

// -----------------------------------------------------------------------------

void ObsWindTLAD::print(std::ostream & os) const {
  os << "QG wind components (u and v) observation operator TL/AD";
}

// -----------------------------------------------------------------------------

}  // namespace qg
