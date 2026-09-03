/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

namespace atlas {
class FunctionSpace;
}  // namespace atlas

namespace eckit {
namespace mpi {
class Comm;
}  // namespace mpi
}  // namespace eckit

namespace oops {

class ProximitySearch;

// -----------------------------------------------------------------------------

/// \brief Partition an interpolation target (lat, lon) to the MPI task owning the nearest
///        source-grid point. Implemented as a handle to a shared ProximitySearch.
class SourceProximityPartitioner {
 public:
  explicit SourceProximityPartitioner(const ProximitySearch & search) : search_(&search) {}
  int interpolatingTask(double lat, double lon) const;

 private:
  const ProximitySearch * search_;
};

/// \brief Builds a SourceProximityPartitioner for the given (\p fspace, \p comm).
///        The first build constructs the underlying cached object, which is a collective
///        operation, so call this uniformly across ranks.
SourceProximityPartitioner makeSourceProximityPartitioner(const atlas::FunctionSpace & fspace,
                                                          const eckit::mpi::Comm & comm);

// -----------------------------------------------------------------------------

}  // namespace oops
