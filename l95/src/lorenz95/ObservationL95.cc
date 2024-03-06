/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/ObservationL95.h"

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "lorenz95/GomL95.h"
#include "lorenz95/L95TraitsFwd.h"
#include "lorenz95/LocsL95.h"
#include "lorenz95/ObsBias.h"
#include "lorenz95/ObsDiags1D.h"
#include "lorenz95/ObsVec1D.h"
#include "oops/base/Locations.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------

ObservationL95::ObservationL95(const ObsTable & ot, const eckit::Configuration &)
  : obsdb_(ot), inputs_(std::vector<std::string>{"x"})
{}

// -----------------------------------------------------------------------------

ObservationL95::~ObservationL95() {}

// -----------------------------------------------------------------------------

void ObservationL95::simulateObs(const GomL95 & gom, ObsVec1D & ovec,
                                 const ObsBias & bias,
                                 const QCFlags_ & qc_flags,
                                 ObsVec1D &, ObsDiags1D &) const {
  for (size_t jj = 0; jj < gom.size(); ++jj) {
    ovec[jj] = gom[jj] + bias.value();
  }
}

// -----------------------------------------------------------------------------

oops::Locations<L95ObsTraits> ObservationL95::locations() const {
  oops::SampledLocations<L95ObsTraits> sampledLocations(
        std::make_unique<LocsL95>(obsdb_.locations(), obsdb_.times()));
  return oops::Locations<L95ObsTraits>(std::move(sampledLocations));
}

// -----------------------------------------------------------------------------

void ObservationL95::print(std::ostream & os) const {
  os << "Lorenz 95: Identity obs operator";
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
