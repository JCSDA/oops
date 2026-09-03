/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/generic/ProximitySearch.h"

#include <cmath>
#include <memory>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

#include "atlas/array.h"
#include "atlas/grid/Grid.h"
#include "atlas/util/Point.h"

#include "eckit/mpi/Comm.h"

#include "oops/mpi/mpi.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------

ProximitySearch::ProximitySearch(const atlas::FunctionSpace & fspace,
                                 const eckit::mpi::Comm & comm):
  comm_(comm), earth_(atlas::util::Earth::radius()), globalNodeTree_(earth_)
{
  // Set default communicator name
  eckit::mpi::setCommDefault(comm_.name().c_str());

  // Exit early if receiving an uninitialized or mesh-less FunctionSpace, because the
  // proximity cache can't be built in this case.
  if (!fspace) {
    Log::info() << "ProximitySearch received uninitialized FunctionSpace,"
      << " so skipping set up of the proximity cache." << std::endl;
    return;
  } else if (fspace.type() == "PointCloud" || fspace.type() == "Spectral") {
    Log::info() << "ProximitySearch received FunctionSpace " << fspace.type()
      << ", so skipping set up of the proximity cache." << std::endl;
    return;
  } else {
    // Check for custom MPI partitions where some tasks handle zero points
    // This configuration could likely be supported after adding new logic to skip/handle work on
    // the zero-size meshes that arise, but we haven't done that yet as it's a bit of an edge case.
    int min_nb_cells = fspace.size();
    comm.allReduceInPlace(min_nb_cells, eckit::mpi::min());
    if (min_nb_cells == 0) {
      Log::info() << "ProximitySearch received FunctionSpace " << fspace.type()
        << " using a distribution where some MPI tasks own zero points"
        << ", so skipping set up of the proximity cache." << std::endl;
      return;
    }
  }

  setGlobalTree(fspace);
}

// -----------------------------------------------------------------------------

int ProximitySearch::taskOwningClosestPoint(const double lat, const double lon) const {
  ASSERT(!globalNodeTree_.empty());
  atlas::PointLonLat target(lon, lat);
  target.normalise();
  const int itask = globalNodeTree_.closestPoint(target).payload().first;
  ASSERT(itask >= 0 && (size_t)itask < comm_.size());
  return itask;
}

// -----------------------------------------------------------------------------

std::optional<int> ProximitySearch::closestPointWithinRadius(const double lat,
    const double lon, const double radius) const {
  ASSERT(!globalNodeTree_.empty());
  ASSERT(radius > 0.0);
  atlas::PointLonLat target(lon, lat);
  target.normalise();
  const auto points = globalNodeTree_.closestPointsWithinRadius(target, radius);
  if (points.size() > 0) {
    ASSERT(points[0].payload().first == static_cast<atlas::idx_t>(comm_.rank()));
    return static_cast<int>(points[0].payload().second);
  } else {
    return std::nullopt;
  }
}

// -----------------------------------------------------------------------------

void ProximitySearch::setGlobalTree(const atlas::FunctionSpace & fspace) {
  ASSERT(globalNodeTree_.empty());
  util::Timer timer("oops::ProximitySearch", "setGlobalTree");

  // Count number of owned points
  const auto lonlat_view = atlas::array::make_view<double, 2>(fspace.lonlat());
  const auto ghost_view = atlas::array::make_view<int, 1>(fspace.ghost());
  const size_t nb_owned = [&ghost_view]() {
    int result = 0;
    for (atlas::idx_t jj = 0; jj < ghost_view.shape(0); ++jj) {
      if (ghost_view(jj) == 0) {
        ++result;
      }
    }
    return result;
  }();

  // Copy owned points into local buffer:
  // lon, lat, and task-local index (as double) into FunctionSpace
  std::vector<double> lonlatidx;
  lonlatidx.reserve(3 * nb_owned);
  for (atlas::idx_t jj = 0; jj < ghost_view.shape(0); ++jj) {
    if (ghost_view(jj) == 0) {
      lonlatidx.emplace_back(lonlat_view(jj, 0));
      lonlatidx.emplace_back(lonlat_view(jj, 1));
      lonlatidx.emplace_back(static_cast<double>(jj));
    }
  }
  ASSERT(lonlatidx.size() == 3 * nb_owned);

  // Collect global grid lats, lons, task-local index
  const size_t nb_tasks = comm_.size();
  std::vector<size_t> sizes(nb_tasks);
  comm_.allGather(nb_owned, sizes.begin(), sizes.end());

  size_t nb_global = 0;
  for (size_t jtask = 0; jtask < nb_tasks; ++jtask) {
    nb_global += sizes[jtask];
  }

  std::vector<double> lonlatidx_global(3 * nb_global);
  mpi::allGatherv(comm_, lonlatidx, lonlatidx_global);

  // Arrange coordinates and {task, task-local index} pair for kd-tree
  std::vector<atlas::PointLonLat> nodes;
  std::vector<std::pair<atlas::idx_t, atlas::idx_t>> payloads;
  nodes.reserve(nb_global);
  payloads.reserve(nb_global);
  int counter = 0;
  for (size_t jtask = 0; jtask < nb_tasks; ++jtask) {
    for (size_t jj = 0; jj < sizes[jtask]; ++jj) {
      nodes.emplace_back(lonlatidx_global[3 * counter], lonlatidx_global[3 * counter + 1]);
      payloads.emplace_back(static_cast<atlas::idx_t>(jtask),
                            static_cast<atlas::idx_t>(lonlatidx_global[3 * counter + 2]));
      ++counter;
    }
  }
  ASSERT(nodes.size() == nb_global);
  ASSERT(payloads.size() == nb_global);
  ASSERT(counter == nb_global);

  // Hacky step:
  // When the source grid is singular and has multiple coinciding grid points, (e.g., a regular
  // latlon grid at the poles), then many nodes of globalNodeTree_ will be at the same physical
  // location, and the KD-tree no longer becomes useful for proximity-based search. Thus, we need
  // to slightly spread out degenerate grid points. Because this is (so far) a rare scenario in
  // JEDI, and because it's unclear how to correctly handle degenerate points in full generality,
  // we write code below to address specific pathological cases.
  //
  // Pathological case of a structured grid with degenerate points at the poles -- shift points
  // away from the poles by a small amount:
  if (fspace.type() == "StructuredColumns") {
    const atlas::functionspace::StructuredColumns structuredcolumns(fspace);
    const atlas::RegularGrid rg(structuredcolumns.grid());
    if (rg) {
      const double eps_check = 1e-14;
      const double shift = 1e-6;  // large epsilon to be robust to trigonometry
      for (auto & lonlat : nodes) {
        double & lat = lonlat[1];
        if (std::abs(lat - 90.0) < eps_check) {
          lat = (90.0 - shift);
        } else if (std::abs(lat + 90.0) < eps_check) {
          lat = -(90.0 - shift);
        }
      }
    }
  }

  // Create global kd-tree
  globalNodeTree_.build(nodes, payloads);
}

// -----------------------------------------------------------------------------

namespace {

/// Flyweight factory to manage creation and retrieval of ProximitySearch objects
class ProximitySearchFactory {
 public:
  static const ProximitySearch & get(const atlas::FunctionSpace & fspace,
                                     const eckit::mpi::Comm & comm) {
    // The cache's payload is MPI task ownership, so it depends on the source grid AND its
    // distribution; the key includes the communicator name and a lon/lat hash in addition to
    // the grid UID -- matching util::CommRedistributionRepository.
    const std::string sep = "_";
    const std::string key = comm.name() + sep + util::getGridUid(fspace) + sep
                            + util::getLonLatHash(fspace);

    auto & map = getInstance().instances_;
    const auto it = map.find(key);
    if (it != map.end()) {
      return *(it->second);
    }
    const auto inserted = map.emplace(key, std::make_unique<ProximitySearch>(fspace, comm));
    Log::debug() << "ProximitySearchFactory: created new entry for key '" << key << "'."
                 << std::endl;
    return *(inserted.first->second);
  }

 private:
  ProximitySearchFactory() {}
  ProximitySearchFactory(const ProximitySearchFactory &) = delete;
  ProximitySearchFactory & operator=(const ProximitySearchFactory &) = delete;

  static ProximitySearchFactory & getInstance() {
    static ProximitySearchFactory theInstance;
    return theInstance;
  }

  std::unordered_map<std::string, std::unique_ptr<ProximitySearch>> instances_;
};

}  // namespace

const ProximitySearch & getProximitySearch(const atlas::FunctionSpace & fspace,
                                           const eckit::mpi::Comm & comm) {
  return ProximitySearchFactory::get(fspace, comm);
}

// -----------------------------------------------------------------------------

}  // namespace oops
