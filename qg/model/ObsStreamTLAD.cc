/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ObsStreamTLAD.h"

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
static ObsOpTLADMaker<ObsStreamTLAD> makerStreamTL_("Stream");
// -----------------------------------------------------------------------------

ObsStreamTLAD::ObsStreamTLAD(const ObsSpaceQG & odb, const eckit::Configuration & config)
  : obsdb_(odb), varin_(std::vector<std::string>{"x", "z"})
{
  oops::Log::trace() << "ObsStreamTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsStreamTLAD::setTrajectory(const GomQG &, const ObsBias &, const QCFlags_ &) {}

// -----------------------------------------------------------------------------

void ObsStreamTLAD::simulateObsTL(const GomQG & gom, ObsVecQG & ovec,
                                  const ObsBiasIncrement & bias, const QCFlags_ & qc_flags) const {
  qg_stream_equiv_tl_f90(obsdb_.toFortran(), gom.toFortran(), ovec.toFortran(), bias.stream());
}

// -----------------------------------------------------------------------------

void ObsStreamTLAD::simulateObsAD(GomQG & gom, const ObsVecQG & ovec,
                                  ObsBiasIncrement & bias, const QCFlags_ & qc_flags) const {
  qg_stream_equiv_ad_f90(obsdb_.toFortran(), gom.toFortran(), ovec.toFortran(), bias.stream());
}

// -----------------------------------------------------------------------------

void ObsStreamTLAD::print(std::ostream & os) const {
  os << "QG Stream observation operator TL/AD";
}

// -----------------------------------------------------------------------------

}  // namespace qg
