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

#include <vector>

#include "eckit/config/Configuration.h"
#include "model/GomQG.h"
#include "model/ObsBias.h"
#include "model/ObsSpaceQG.h"
#include "model/ObsVecQG.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
static ObsOpMaker<ObsStreamQG> makerStream_("Stream");
// -----------------------------------------------------------------------------

ObsStreamQG::ObsStreamQG(const ObsSpaceQG & odb, const eckit::Configuration & config)
  : keyOperStrm_(0), obsdb_(odb), varin_(std::vector<std::string>{"x"})
{
  qg_stream_setup_f90(keyOperStrm_, config);
  oops::Log::trace() << "ObsStreamQG created." << std::endl;
}

// -----------------------------------------------------------------------------

ObsStreamQG::~ObsStreamQG() {
  qg_stream_delete_f90(keyOperStrm_);
  oops::Log::trace() << "ObsStreamQG destructed." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsStreamQG::simulateObs(const GomQG & gom, ObsVecQG & ovec,
                              const ObsBias & bias) const {
  qg_stream_equiv_f90(gom.toFortran(), ovec.toFortran(), bias.stream());
}

// -----------------------------------------------------------------------------

LocationsQG * ObsStreamQG::locations(const util::DateTime & t1, const util::DateTime & t2) const {
  return obsdb_.locations(t1, t2);
}

// -----------------------------------------------------------------------------

void ObsStreamQG::print(std::ostream & os) const {
  os << "ObsStreamQG::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace qg
