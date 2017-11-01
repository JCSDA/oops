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

#include <string>
#include <vector>

#include "util/Logger.h"
#include "eckit/config/Configuration.h"
#include "util/DateTime.h"

#include "lorenz95/GomL95.h"
#include "lorenz95/NoVariables.h"
#include "lorenz95/ObsBias.h"
#include "lorenz95/ObsVec1D.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------

ObservationL95::ObservationL95(ObsTable & ot, const eckit::Configuration & conf)
  : obsdb_(ot), inputs_(new NoVariables())
{}

// -----------------------------------------------------------------------------

ObservationL95::~ObservationL95() {}

// -----------------------------------------------------------------------------

void ObservationL95::obsEquiv(const GomL95 & gom, ObsVec1D & ovec,
                              const ObsBias & bias) const {
  for (int jj = 0; jj < gom.nobs(); ++jj) {
    const int ii = gom.getindx(jj);
    ovec(ii)=gom[jj] + bias.value();
  }
}

// -----------------------------------------------------------------------------

void ObservationL95::generateObsError(const eckit::Configuration & conf) {
  const double err = conf.getDouble("obs_error");
  std::vector<double> obserr(obsdb_.nobs());
  for (unsigned int jj = 0; jj < obserr.size(); ++jj) {
    obserr[jj] = err;
  }
  obsdb_.putdb("ObsErr", obserr);
}

// -----------------------------------------------------------------------------

void ObservationL95::print(std::ostream & os) const {
  os << "ObservationL95::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
