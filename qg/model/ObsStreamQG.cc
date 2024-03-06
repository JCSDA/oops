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
#include "model/LocationsQG.h"
#include "model/ObsBias.h"
#include "model/ObsSpaceQG.h"
#include "model/ObsVecQG.h"
#include "model/QgTraitsFwd.h"
#include "oops/base/Locations.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
static ObsOpMaker<ObsStreamQG> makerStream_("Stream");
// -----------------------------------------------------------------------------

ObsStreamQG::ObsStreamQG(const ObsSpaceQG & odb, const eckit::Configuration & config)
  : obsdb_(odb), varin_(std::vector<std::string>{"x", "z"})
{
  oops::Log::trace() << "ObsStreamQG created." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsStreamQG::simulateObs(const GomQG & gom, ObsVecQG & ovec,
                              const ObsBias & bias,
                              const QCFlags_ & qc_flags) const {
  qg_stream_equiv_f90(obsdb_.toFortran(), gom.toFortran(), ovec.toFortran(), bias.stream());
}

// -----------------------------------------------------------------------------

ObsOpBaseQG::Locations_ ObsStreamQG::locations() const {
  typedef oops::SampledLocations<QgObsTraits> SampledLocations_;
  return Locations_(SampledLocations_(obsdb_.locations()));
}

// -----------------------------------------------------------------------------

void ObsStreamQG::print(std::ostream & os) const {
  os << "QG Stream observation operator TL/AD";
}

// -----------------------------------------------------------------------------

}  // namespace qg
