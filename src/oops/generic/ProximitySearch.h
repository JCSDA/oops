/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <optional>
#include <utility>

#include "atlas/functionspace.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"

namespace eckit {
namespace mpi {
class Comm;
}  // namespace mpi
}  // namespace eckit

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Proximity-based search over a grid: given an target point (lat, lon), find which
///        MPI task owns the nearest source point, or the local index of the nearest source
///        point within a radius.
///
/// Construction is MPI-collective (allGatherv), so it must be reached uniformly across all
/// ranks of \p comm.
class ProximitySearch {
 public:
  ProximitySearch(const atlas::FunctionSpace &, const eckit::mpi::Comm &);

  ~ProximitySearch() = default;

  ProximitySearch(const ProximitySearch &) = delete;
  ProximitySearch & operator=(const ProximitySearch &) = delete;

  /// Returns the MPI task that owns the source point globally nearest to (lat, lon).
  int taskOwningClosestPoint(const double lat, const double lon) const;

  /// Returns the task-local index into the FunctionSpace of the globally-nearest
  /// grid point to (lat, lon). Aborts if the cache is empty or if the nearest point
  /// is not owned by this MPI rank. It is assumed that the taskOwningClosestPoint method will
  /// have been previously used, so that the call to this method is placed on the
  /// correct MPI task.
  /// Note that `radius` is a chord distance in meters.
  std::optional<int> closestPointWithinRadius(double lat, double lon, double radius) const;

 private:
  void setGlobalTree(const atlas::FunctionSpace &);

  const eckit::mpi::Comm & comm_;
  const atlas::Geometry earth_;
  // JEDI grid nodes = model cell-centers; payload = {owning task, task-local index}
  atlas::util::KDTree<std::pair<atlas::idx_t, atlas::idx_t>> globalNodeTree_;
};

// -----------------------------------------------------------------------------

/// \brief Returns the ProximitySearch for the source geometry described by
///        (\p fspace, \p comm), building and caching it on first request.
///
/// The cache is held in a process-local registry keyed by grid identity (communicator name
/// + grid UID + lon/lat hash), so repeated interpolators/orchestrators on the same source
/// geometry share a single build. Because construction is MPI-collective, this must be called
/// uniformly across all ranks of \p comm.
///
/// Computing the registry's key can be up to O(grid size), so this function should be
/// called outside of any loops, and the returned reference can be queried within a loop.
const ProximitySearch & getProximitySearch(const atlas::FunctionSpace & fspace,
                                           const eckit::mpi::Comm & comm);

// -----------------------------------------------------------------------------

}  // namespace oops
