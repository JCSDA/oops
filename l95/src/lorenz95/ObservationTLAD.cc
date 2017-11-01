/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/ObservationTLAD.h"

#include <string>

#include "eckit/config/Configuration.h"
#include "util/Logger.h"
#include "util/DateTime.h"

#include "lorenz95/GomL95.h"
#include "lorenz95/NoVariables.h"
#include "lorenz95/ObsBiasCorrection.h"
#include "lorenz95/ObsTable.h"
#include "lorenz95/ObsVec1D.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------

ObservationTLAD::ObservationTLAD(const ObsTable &)
  : inputs_(new NoVariables())
{}

// -----------------------------------------------------------------------------

ObservationTLAD::~ObservationTLAD() {}

// -----------------------------------------------------------------------------

void ObservationTLAD::setTrajectory(const GomL95 &, const ObsBias &) {}

// -----------------------------------------------------------------------------

void ObservationTLAD::obsEquivTL(const GomL95 & gom, ObsVec1D & ovec,
                                 const ObsBiasCorrection & bias) const {
  for (int jj = 0; jj < gom.nobs(); ++jj) {
    const int ii = gom.getindx(jj);
    ovec(ii) = gom[jj] + bias.value();
  }
}

// -----------------------------------------------------------------------------

void ObservationTLAD::obsEquivAD(GomL95 & gom, const ObsVec1D & ovec,
                                 ObsBiasCorrection & bias) const {
  for (int jj = 0; jj < gom.nobs(); ++jj) {
    const int ii = gom.getindx(jj);
    gom[jj] = ovec(ii);
    bias.value() += ovec(ii);
  }
}

// -----------------------------------------------------------------------------

void ObservationTLAD::print(std::ostream & os) const {
  os << "ObservationTLAD::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
