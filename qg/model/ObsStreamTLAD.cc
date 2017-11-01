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

#include "eckit/config/Configuration.h"
#include "model/GomQG.h"
#include "model/ObsBias.h"
#include "model/ObsBiasIncrement.h"
#include "model/ObsSpaceQG.h"
#include "model/ObsVecQG.h"
#include "model/QgFortran.h"
#include "model/VariablesQG.h"
#include "util/Logger.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
static oops::LinearObsOpMaker<QgTraits, ObsStreamTLAD> makerStreamTL_("Stream");
// -----------------------------------------------------------------------------

ObsStreamTLAD::ObsStreamTLAD(const ObsSpaceQG &, const eckit::Configuration & config)
  : keyOperStrm_(0), varin_()
{
  const eckit::Configuration * configc = &config;
  qg_stream_setup_f90(keyOperStrm_, &configc);
  int keyVarin;
  qg_obsoper_inputs_f90(keyOperStrm_, keyVarin);
  varin_.reset(new VariablesQG(keyVarin));
  oops::Log::trace() << "ObsStreamTLAD created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsStreamTLAD::~ObsStreamTLAD() {
  oops::Log::trace() << "ObsStreamTLAD destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsStreamTLAD::setTrajectory(const GomQG &, const ObsBias &) {}

// -----------------------------------------------------------------------------

void ObsStreamTLAD::obsEquivTL(const GomQG & gom, ObsVecQG & ovec,
                               const ObsBiasIncrement & bias) const {
  qg_stream_equiv_tl_f90(gom.toFortran(), ovec.toFortran(), bias.stream());
}

// -----------------------------------------------------------------------------

void ObsStreamTLAD::obsEquivAD(GomQG & gom, const ObsVecQG & ovec,
                               ObsBiasIncrement & bias) const {
  qg_stream_equiv_ad_f90(gom.toFortran(), ovec.toFortran(), bias.stream());
}

// -----------------------------------------------------------------------------

void ObsStreamTLAD::print(std::ostream & os) const {
  os << "ObsStreamTLAD::print not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace qg
