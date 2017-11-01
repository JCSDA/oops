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
#include "lorenz95/ObsTable.h"
#include "util/DateTime.h"

namespace lorenz95 {

// -----------------------------------------------------------------------------

LocsL95::LocsL95(const ObsTable & ot,
                 const util::DateTime & t1, const util::DateTime & t2) {
  locs_ = ot.locations(t1, t2);
}

// -----------------------------------------------------------------------------

void LocsL95::print(std::ostream & os) const {
  os << "LocsL95::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
