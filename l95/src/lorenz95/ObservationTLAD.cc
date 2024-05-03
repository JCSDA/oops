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
#include <vector>

#include "lorenz95/GomL95.h"
#include "lorenz95/ObsBiasCorrection.h"
#include "lorenz95/ObsVec1D.h"
#include "oops/base/Variables.h"
#include "oops/util/missingValues.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------

ObservationTLAD::ObservationTLAD(const ObsTable &, const eckit::Configuration &)
  : inputs_(std::vector<std::string>{"x"})
{}

// -----------------------------------------------------------------------------

void ObservationTLAD::setTrajectory(const GomL95 &, const ObsBias &, const QCFlags_ &) {}

// -----------------------------------------------------------------------------

void ObservationTLAD::simulateObsTL(const GomL95 & gom, ObsVec1D & ovec,
                                    const ObsBiasCorrection & bias,
                                    const QCFlags_ & qc_flags) const {
  for (size_t jj = 0; jj < gom.size(); ++jj) {
    ovec[jj] = gom[jj] + bias.value();
  }
}

// -----------------------------------------------------------------------------

void ObservationTLAD::simulateObsAD(GomL95 & gom, const ObsVec1D & ovec,
                                    ObsBiasCorrection & bias,
                                    const QCFlags_ & qc_flags) const {
  const double missing = util::missingValue<double>();
  for (size_t jj = 0; jj < gom.size(); ++jj) {
    if (ovec[jj] != missing) {
      gom[jj] += ovec[jj];
      bias.value() += ovec[jj];
    }
  }
}

// -----------------------------------------------------------------------------

void ObservationTLAD::print(std::ostream & os) const {
  os << "Lorenz 95: Identity linear obs operator";
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
