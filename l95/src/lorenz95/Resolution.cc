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
#include <vector>

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------
std::vector<int> Resolution::getDims() const {
  std::vector<int> dims(1);
  dims[0] = resol_;
  return dims;
}
// -----------------------------------------------------------------------------
std::vector<double> Resolution::getLats() const {
  std::vector<double> lats(resol_);
  for (int jj = 0; jj < resol_; ++jj) lats[jj] = 0.0;
  return lats;
}
// -----------------------------------------------------------------------------
std::vector<double> Resolution::getLons() const {
  std::vector<double> lons(resol_);
  double dx = 360.0 / resol_;
  for (int jj = 0; jj < resol_; ++jj) lons[jj] = dx * jj;
  return lons;
}
// -----------------------------------------------------------------------------
std::vector<double> Resolution::getLevs() const {
  std::vector<double> levs(1);
  levs[0] = 0.0;
  return levs;
}
// -----------------------------------------------------------------------------
std::vector<double> Resolution::getArea() const {
  std::vector<double> area(1);
  area[0] = 1.0;
  return area;
}
// -----------------------------------------------------------------------------
std::vector<int> Resolution::getMask(const int &) const {
  std::vector<int> mask(resol_);
  for (int jj = 0; jj < resol_; ++jj) mask[jj] = 1;
  return mask;
}
// -----------------------------------------------------------------------------

}  // namespace lorenz95
