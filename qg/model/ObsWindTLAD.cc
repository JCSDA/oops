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

ObsWindTLAD::ObsWindTLAD(const ObsSpaceQG & odb, const int & keyOperWind)
  : keyOperWind_(keyOperWind), varin_()
{
  int keyVarin;
  qg_obsoper_inputs_f90(keyOperWind_, keyVarin);
  varin_.reset(new VariablesQG(keyVarin));
  Log::trace() << "ObsWindTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsWindTLAD::~ObsWindTLAD() {
  Log::trace() << "ObsWindTLAD destrcuted" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsWindTLAD::setTrajectory(const GomQG &, const ObsBias &) {}

// -----------------------------------------------------------------------------

void ObsWindTLAD::obsEquivTL(const GomQG & gom, ObsVecQG & ovec,
                             const ObsBiasIncrement & bias) const {
  qg_wind_equiv_tl_f90(gom.toFortran(), ovec.toFortran(), bias.wind());
}

// -----------------------------------------------------------------------------

void ObsWindTLAD::obsEquivAD(GomQG & gom, const ObsVecQG & ovec,
                             ObsBiasIncrement & bias) const {
  qg_wind_equiv_ad_f90(gom.toFortran(), ovec.toFortran(), bias.wind());
}

// -----------------------------------------------------------------------------

}  // namespace qg
