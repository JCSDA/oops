/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/GomL95.h"

#include <cstdlib>
#include <vector>

#include "lorenz95/ObsTable.h"
#include "lorenz95/NoVariables.h"
#include "util/Logger.h"

namespace lorenz95 {
  class Resolution;

// -----------------------------------------------------------------------------
GomL95::GomL95(const ObsTable & ot, const NoVariables &,
               const util::DateTime & t1, const util::DateTime & t2,
               const Resolution &)
  : size_(0), iobs_(), locval_(), current_(0)
{
  iobs_ = ot.timeSelect(t1, t2);
  size_ = iobs_.size();
  locval_.resize(size_);
}
// -----------------------------------------------------------------------------
GomL95::~GomL95() {}
// -----------------------------------------------------------------------------
void GomL95::zero() {
  for (int jj = 0; jj < size_; ++jj) locval_[jj] = 0.0;
}
// -----------------------------------------------------------------------------
double GomL95::dot_product_with(const GomL95 & gom) const {
  double zz = 0.0;
  for (int jj = 0; jj < size_; ++jj) zz += locval_[jj] * gom.locval_[jj];
  return zz;
}
// -----------------------------------------------------------------------------
void GomL95::print(std::ostream & os) const {
  os << "GomL95::print not implemented";
}
// -----------------------------------------------------------------------------

}  // namespace lorenz95
