/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "model/LocationsQG.h"
#include "model/QgFortran.h"

namespace qg {

// -------------------------------------------------------------------------
LocationsQG::LocationsQG(const eckit::Configuration & config, const eckit::mpi::Comm &) {
  qg_locs_create_f90(keyLocs_);
  if (config.has("lats") || config.has("Nrandom")) {
    std::vector<double> lons = config.getDoubleVector("lons");
    std::vector<double> lats = config.getDoubleVector("lats");
    std::vector<double> z  = config.getDoubleVector("z");

    ASSERT(lons.size() == lats.size());
    ASSERT(lons.size() == z.size());
    const unsigned int nlocs = lons.size();

    qg_locs_test_f90(keyLocs_, config, nlocs, &lons[0], &lats[0], &z[0]);
  }
}
// -------------------------------------------------------------------------
int LocationsQG::size() const {
  int nobs = 0;
  qg_locs_nobs_f90(keyLocs_, nobs);
  return nobs;
}
// -------------------------------------------------------------------------
void LocationsQG::print(std::ostream & os) const {
  int nobs = 0;
  qg_locs_nobs_f90(keyLocs_, nobs);
  double lon = 0.0;
  double lat = 0.0;
  double z = 0.0;
  for (size_t jj=0; jj < static_cast<size_t>(nobs); ++jj) {
    qg_locs_element_f90(keyLocs_, jj, lon, lat, z);
    os << "location " << jj << std::setprecision(2) << ": lon = " << lon
       << ", lat = " << lat << ", z = " << z << std::endl;
  }
}
// -------------------------------------------------------------------------
}  // namespace qg
