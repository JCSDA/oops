/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/LocsL95.h"

#include "eckit/config/LocalConfiguration.h"
#include "lorenz95/ObsTable.h"
#include "util/DateTime.h"
#include "util/Logger.h"

namespace lorenz95 {

// -----------------------------------------------------------------------------

LocsL95::LocsL95(const std::vector<int> & indx, const std::vector<double> & locs)
 : indx_(indx), locs_(locs)
{
  ASSERT(indx_.size() == locs_.size());
}

// -----------------------------------------------------------------------------

LocsL95::LocsL95(const eckit::Configuration & conf) : indx_(), locs_() {
  const double zz = conf.getDouble("position");
  ASSERT(zz >= 0.0 && zz <= 1.0);
  locs_.push_back(zz);
  indx_.push_back(1);
  oops::Log::trace() << "LocsL95::LocsL95 created" << std::endl;
}

// -----------------------------------------------------------------------------

void LocsL95::print(std::ostream & os) const {
  os << locs_.size();
  if (locs_.size() > 0) os << " " << locs_.at(1);
  if (locs_.size() > 1) os << " " << locs_.at(locs_.size() - 1);
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
