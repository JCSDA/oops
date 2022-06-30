/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2021-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/LocsL95.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "lorenz95/ObsTable.h"
#include "oops/util/DateTime.h"

namespace lorenz95 {

// -----------------------------------------------------------------------------

LocsL95::LocsL95(const std::vector<double> & locs, const std::vector<util::DateTime> & times)
  : locs_(locs), dummy_(locs.size(), 0.0), times_(times)
{
  ASSERT(locs_.size() == times_.size());
}

// -----------------------------------------------------------------------------

LocsL95::LocsL95(const eckit::Configuration & conf, const eckit::mpi::Comm &)
  : locs_(), dummy_(), times_() {
  conf.get("positions", locs_);
  const util::DateTime time(conf.getString("time"));
  for (size_t jj = 0; jj < locs_.size(); ++jj) {
    ASSERT(locs_.at(jj) >= 0.0 && locs_.at(jj) <= 1.0);
    times_.push_back(time);
    dummy_.push_back(0.0);
  }
}

// -----------------------------------------------------------------------------

void LocsL95::print(std::ostream & os) const {
  os << locs_.size();
  if (locs_.size() > 0) os << " " << locs_.at(0);
  if (locs_.size() > 1) os << " " << locs_.at(locs_.size() - 1);
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
