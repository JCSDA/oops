/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/Resolution.h"

#include <string>
#include <vector>

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------
Iterator Resolution::begin() const {
  return Iterator(*this, 0);
}
// -----------------------------------------------------------------------------
Iterator Resolution::end() const {
  return Iterator(*this, resol_);
}
// -----------------------------------------------------------------------------
std::vector<double> Resolution::verticalCoord(std::string & vcUnits) const {
  std::vector<double> vc(1, 1.0);
  return vc;
}
// -----------------------------------------------------------------------------
std::vector<size_t> Resolution::variableSizes(const oops::Variables & vars) const {
  std::vector<size_t> sizes(vars.size(), 1);
  return sizes;
}
// -----------------------------------------------------------------------------
void Resolution::latlon(std::vector<double> & lats, std::vector<double> & lons, const bool) const {
  const double dx = 1.0 / static_cast<double>(resol_);
  lats.resize(resol_);
  lons.resize(resol_);
  for (size_t jj = 0; jj < (size_t)resol_; ++jj) {
    lons[jj] = static_cast<double>(jj) * dx;
    lats[jj] = 0.0;
  }
}
// -----------------------------------------------------------------------------

}  // namespace lorenz95
