/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ObsWSpeedQG.h"

#include <vector>

#include "eckit/config/Configuration.h"
#include "model/GomQG.h"
#include "model/ObsBias.h"
#include "model/ObsSpaceQG.h"
#include "model/ObsVecQG.h"
#include "model/QgFortran.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------
static ObsOpMaker<ObsWSpeedQG> makerWSpeed_("WSpeed");
// -----------------------------------------------------------------------------

ObsWSpeedQG::ObsWSpeedQG(const ObsSpaceQG & odb, const eckit::Configuration & config)
  : obsdb_(odb), varin_(std::vector<std::string>{"u", "v"})
{
  oops::Log::trace() << "ObsWSpeedQG created." << std::endl;
}

// -----------------------------------------------------------------------------

void ObsWSpeedQG::simulateObs(const GomQG & gom, ObsVecQG & ovec,
                              const ObsBias & bias) const {
  qg_wspeed_equiv_f90(gom.toFortran(), ovec.toFortran(), bias.wspd());
}

// -----------------------------------------------------------------------------

std::unique_ptr<LocationsQG> ObsWSpeedQG::locations() const {
  return obsdb_.locations();
}

// -----------------------------------------------------------------------------

void ObsWSpeedQG::print(std::ostream & os) const {
  os << "QG wind speed observation operator TL/AD";
}

// -----------------------------------------------------------------------------

}  // namespace qg
