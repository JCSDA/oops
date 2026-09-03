/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/generic/SourceProximityPartitioner.h"

#include "oops/generic/ProximitySearch.h"

namespace oops {

// -----------------------------------------------------------------------------

int SourceProximityPartitioner::interpolatingTask(double lat, double lon) const {
  return search_->taskOwningClosestPoint(lat, lon);
}

// -----------------------------------------------------------------------------

SourceProximityPartitioner makeSourceProximityPartitioner(const atlas::FunctionSpace & fspace,
                                                          const eckit::mpi::Comm & comm) {
  return SourceProximityPartitioner{getProximitySearch(fspace, comm)};
}

// -----------------------------------------------------------------------------

}  // namespace oops
