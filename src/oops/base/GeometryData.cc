/*
 * (C) Copyright 2022- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/GeometryData.h"

#include "oops/external/stripack/stripack.h"
#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------

GeometryData::GeometryData(const atlas::FunctionSpace & fspace, const atlas::FieldSet & fset,
                           const bool topdown, const eckit::mpi::Comm & comm):
  fspace_(fspace), fset_(fset), comm_(&comm), topdown_(topdown),
  earth_(atlas::util::Earth::radius()), localTree_(earth_), globalTree_(earth_),
  loctree_(false), glotree_(false),
  unitsphere_(1.0), triangulation_(nullptr)
{}

// -----------------------------------------------------------------------------

GeometryData::~GeometryData() = default;

// -----------------------------------------------------------------------------

int GeometryData::closestTask(const double lat, const double lon) const {
  ASSERT(glotree_);
  atlas::PointLonLat obsloc(lon, lat);
  obsloc.normalise();
  const int itask = globalTree_.closestPoint(obsloc).payload();
  ASSERT(itask >= 0 && (size_t)itask < comm_->size());
  return itask;
}

// -----------------------------------------------------------------------------

atlas::util::KDTree<size_t>::ValueList GeometryData::closestPoints(const double lat,
                                                                   const double lon,
                                                                   const int npoints) const {
  ASSERT(loctree_);
  atlas::PointLonLat obsloc(lon, lat);
  obsloc.normalise();
  atlas::util::KDTree<size_t>::ValueList neighbours = localTree_.closestPoints(obsloc, npoints);
  return neighbours;
}

// -----------------------------------------------------------------------------

bool GeometryData::containingTriangleAndBarycentricCoords(const double lat, const double lon,
    std::array<int, 3> & indices, std::array<double, 3> & baryCoords) const {
  ASSERT(loctree_);

  // Initialize Triangulation on first use
  if (triangulation_ == nullptr) {
    ASSERT(lats_.size() > 0);  // check temporary coord buffers were initialized
    ASSERT(lats_.size() == lons_.size());
    triangulation_.reset(new stripack::Triangulation(lats_, lons_));
  }

  atlas::PointLonLat ptll(lon, lat);
  ptll.normalise();
  atlas::Point3 pt3;
  unitsphere_.lonlat2xyz(ptll, pt3);
  const std::array<double, 3> coords = {{pt3[0], pt3[1], pt3[2]}};

  // Index of closest neighbor in source-point list
  const int guessIndex = localTree_.closestPoint(ptll).payload();

  const bool validTriangle = triangulation_->containingTriangleAndBarycentricCoords(
      coords, guessIndex, indices, baryCoords);
  return validTriangle;
}

// -----------------------------------------------------------------------------

// Local tree requires lats and lons with halo
void GeometryData::setLocalTree(const std::vector<double> & lats,
                                const std::vector<double> & lons) {
  util::Timer timer("oops::GeometryData", "setLocalTree");
  const size_t npoints = lats.size();
  std::vector<size_t> indx(npoints);
  for (size_t jj = 0; jj < npoints; ++jj) indx[jj] = jj;
  localTree_.build(lons, lats, indx);
  loctree_ = true;

  // Save lats,lons for constructing the triangulation
  // (triangulation is not valid for all geometries; only construct it on demand)
  lats_ = lats;
  lons_ = lons;
}

// -----------------------------------------------------------------------------

// Global tree requires lats and lons without halo
void GeometryData::setGlobalTree(const std::vector<double> & lats,
                                 const std::vector<double> & lons) {
  util::Timer timer("oops::GeometryData", "setGlobalTree");
  const size_t ntasks = comm_->size();

// Local latitudes and longitudes
  const size_t sizel = lats.size();
  std::vector<double> latlon(2*sizel);
  for (size_t jj = 0; jj < sizel; ++jj) {
    latlon[2*jj]   = lats[jj];
    latlon[2*jj+1] = lons[jj];
  }

// Collect global grid lats and lons
  std::vector<size_t> sizes(ntasks);
  comm_->allGather(sizel, sizes.begin(), sizes.end());

  size_t sizeg = 0;
  for (size_t jtask = 0; jtask < ntasks; ++jtask) sizeg += sizes[jtask];

  std::vector<double> latlon_global(2*sizeg);
  mpi::allGatherv(*comm_, latlon, latlon_global);

// Arrange coordinates and task index for kd-tree
  std::vector<double> latglo(sizeg), longlo(sizeg);
  std::vector<int> taskindx(sizeg);
  size_t jglo = 0;
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    for (size_t jj = 0; jj < sizes[jtask]; ++jj) {
      latglo[jglo] = latlon_global[2*jglo];
      longlo[jglo] = latlon_global[2*jglo+1];
      taskindx[jglo] = jtask;
      ++jglo;
    }
  }
  ASSERT(jglo == sizeg);

// Create global kd-tree
  globalTree_.build(longlo, latglo, taskindx);
  glotree_ = true;

  Log::info() << "GeometryData: Global tree size = " << globalTree_.size()
              << ", footprint = " << globalTree_.footprint() << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops
