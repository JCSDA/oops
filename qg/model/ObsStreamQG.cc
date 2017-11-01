/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ObsStreamQG.h"

#include "util/Logger.h"
#include "model/GomQG.h"
#include "model/ObsBias.h"
#include "model/ObsSpaceQG.h"
#include "model/ObsVecQG.h"
#include "model/VariablesQG.h"
#include "eckit/config/Configuration.h"

using oops::Log;

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------

ObsStreamQG::ObsStreamQG(ObsSpaceQG & odb, const eckit::Configuration & config)
  : obsdb_(odb), obsname_("Stream"), varin_()
{
  const eckit::Configuration * configc = &config;
  qg_stream_setup_f90(keyOperStrm_, &configc);
  int keyVarin;
  qg_obsoper_inputs_f90(keyOperStrm_, keyVarin);
  varin_.reset(new VariablesQG(keyVarin));
  Log::trace() << "ObsStreamQG created " << obsname_ << std::endl;
}

// -----------------------------------------------------------------------------

ObsStreamQG::~ObsStreamQG() {
  qg_stream_delete_f90(keyOperStrm_);
}

// -----------------------------------------------------------------------------

void ObsStreamQG::obsEquiv(const GomQG & gom, ObsVecQG & ovec,
                           const ObsBias & bias) const {
  qg_stream_equiv_f90(gom.toFortran(), ovec.toFortran(), bias.stream());
}

// -----------------------------------------------------------------------------

void ObsStreamQG::generateObsError(const eckit::Configuration & conf) {
  const double err = conf.getDouble("obs_error");
  qg_obsdb_seterr_f90(obsdb_.toFortran(), keyOperStrm_, err);
}

// -----------------------------------------------------------------------------

void ObsStreamQG::print(std::ostream & os) const {
  os << "ObsStreamQG::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace qg
