/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2021- UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_GEOMETRY_H_
#define OOPS_BASE_GEOMETRY_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"

#include "oops/interface/Geometry.h"
#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"
#include "oops/util/TypeTraits.h"

namespace eckit {
  class Configuration;
}

namespace oops {

/// \brief Checks whether Geometry has method closestTask. Default: no.
template<class, class = void>
struct HasClosestTask
  : std::false_type {};

/// \brief Checks whether Geometry has method closestTask. Specialization for the case
///        when it does.
template<class Geometry>
struct HasClosestTask<Geometry,
       util::void_t<decltype(std::declval<Geometry>().closestTask(double(), double()))>>
  : std::true_type {};

// -----------------------------------------------------------------------------
/// \brief Geometry class used in oops; subclass of interface class interface::Geometry.
///
/// \details Handles additional MPI communicator parameter in the constructors
/// (for MPI distribution in time, used in oops for 4DEnVar and weak-constraint 4DVar).
/// Adds extra methods that do not need to be implemented in the implementations:
/// - timeComm() (accessor to the MPI communicator in time)
template <typename MODEL>
class Geometry : public interface::Geometry<MODEL> {
  typedef typename MODEL::Geometry              Geometry_;

 public:
  typedef typename interface::Geometry<MODEL>::Parameters_ Parameters_;

  /// Constructor from Parameters and mpi communicators: \p geometry for spatial distribution
  /// (handled by the implementation) and \p time for distribution in time (handled by oops)
  Geometry(const Parameters_ &, const eckit::mpi::Comm & geometry,
           const eckit::mpi::Comm & time = oops::mpi::myself());
  /// Constructor from Configuration and mpi communicators: \p geometry for spatial distribution
  /// (handled by the implementation) and \p time for distribution in time (handled by oops)
  Geometry(const eckit::Configuration &, const eckit::mpi::Comm & geometry,
           const eckit::mpi::Comm & time = oops::mpi::myself());
  /// Constructor from pointer to the MODEL::Geometry (used in 1DVar filter)
  explicit Geometry(std::shared_ptr<const Geometry_>);

  Geometry(const Geometry &) = delete;
  Geometry & operator=(const Geometry &) = delete;

  bool operator==(const Geometry &) const;

  /// Accessor to the MPI communicator for distribution in time
  const eckit::mpi::Comm & timeComm() const {return *timeComm_;}

  /// Returns the MPI task that contains the closest point to the point with
  /// specified \p lat and \p lon.
  ///@{
  /// If MODEL::Geometry has method closestTask implemented, call it.
  template<class Geom = Geometry_>
  typename std::enable_if< HasClosestTask<Geom>::value, int>::type
  closestTask(const double lat, const double lon) const {
    return this->geom_->closestTask(lat, lon);
  }

  /// If MODEL::Geometry doesn't have closestTask implemented,
  /// use a generic implementation.
  template<class Geom = Geometry_>
  typename std::enable_if<!HasClosestTask<Geom>::value, int>::type
  closestTask(const double lat, const double lon) const {
    atlas::PointLonLat obsloc(lon, lat);
    obsloc.normalise();
    const int itask = globalTree_.closestPoint(obsloc).payload();
    ASSERT(itask >= 0 && (size_t)itask < spaceComm_->size());
    return itask;
  }
  ///@}

  atlas::util::KDTree<size_t>::ValueList closestPoints(const double, const double, const int) const;

 private:
  std::unique_ptr<util::Timer> timer_;
  const eckit::mpi::Comm * spaceComm_;  /// pointer to the MPI communicator in space
  const eckit::mpi::Comm * timeComm_;   /// pointer to the MPI communicator in time

  const atlas::Geometry earth_;
  atlas::util::IndexKDTree localTree_;
  atlas::util::IndexKDTree globalTree_;

  void setLocalTree();
  void setGlobalTree();
};

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::Geometry(const eckit::Configuration & config,
                          const eckit::mpi::Comm & geometry, const eckit::mpi::Comm & time):
  interface::Geometry<MODEL>(config, geometry),
  timer_(std::make_unique<util::Timer>("oops::base::Geometry", "Geometry")),
  spaceComm_(&geometry), timeComm_(&time),
  earth_(atlas::util::Earth::radius()), localTree_(earth_), globalTree_(earth_)
{
  this->setLocalTree();
  this->setGlobalTree();
  timer_.reset();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::Geometry(const Parameters_ & parameters,
                          const eckit::mpi::Comm & geometry, const eckit::mpi::Comm & time):
  interface::Geometry<MODEL>(parameters, geometry),
  timer_(std::make_unique<util::Timer>("oops::base::Geometry", "Geometry")),
  spaceComm_(&geometry), timeComm_(&time),
  earth_(atlas::util::Earth::radius()), localTree_(earth_), globalTree_(earth_)
{
  this->setLocalTree();
  this->setGlobalTree();
  timer_.reset();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::Geometry(std::shared_ptr<const Geometry_> ptr):
  interface::Geometry<MODEL>(ptr),
  timer_(std::make_unique<util::Timer>("oops::base::Geometry", "Geometry")),
  spaceComm_(&oops::mpi::world()), timeComm_(&oops::mpi::myself()),
  earth_(atlas::util::Earth::radius()), localTree_(earth_), globalTree_(earth_)
{
  this->setLocalTree();
  this->setGlobalTree();
  timer_.reset();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
bool Geometry<MODEL>::operator==(const Geometry & rhs) const {
  Log::trace() << "Geometry<MODEL>::operator== starting" << std::endl;
  bool eq = (this->geom_ == rhs.geom_ &&
             spaceComm_ == rhs.spaceComm_ && timeComm_ == rhs.timeComm_);
  Log::trace() << "Geometry<MODEL>::operator== done" << std::endl;
  return eq;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
atlas::util::KDTree<size_t>::ValueList
    Geometry<MODEL>::closestPoints(const double lat, const double lon, const int npoints) const {
  atlas::PointLonLat obsloc(lon, lat);
  obsloc.normalise();
  atlas::util::KDTree<size_t>::ValueList neighbours =
                                  localTree_.closestPoints(obsloc, npoints);
  return neighbours;
}

// -----------------------------------------------------------------------------
//  Private methods
// -----------------------------------------------------------------------------

template <typename MODEL>
void Geometry<MODEL>::setLocalTree() {
  std::vector<double> lats;
  std::vector<double> lons;
  this->latlon(lats, lons, true);
  const size_t npoints = lats.size();
  std::vector<size_t> indx(npoints);
  for (size_t jj = 0; jj < npoints; ++jj) indx[jj] = jj;
  localTree_.build(lons, lats, indx);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void Geometry<MODEL>::setGlobalTree() {
  const size_t ntasks = spaceComm_->size();

// Local latitudes and longitudes
  std::vector<double> lats;
  std::vector<double> lons;
  this->latlon(lats, lons, false);
  const size_t sizel = lats.size();
  std::vector<double> latlon(2*sizel);
  for (size_t jj = 0; jj < sizel; ++jj) {
    latlon[2*jj]   = lats[jj];
    latlon[2*jj+1] = lons[jj];
  }

// Collect global grid lats and lons
  std::vector<size_t> sizes(ntasks);
  spaceComm_->allGather(sizel, sizes.begin(), sizes.end());

  size_t sizeg = 0;
  for (size_t jtask = 0; jtask < ntasks; ++jtask) sizeg += sizes[jtask];

  std::vector<double> latlon_global(2*sizeg);
  mpi::allGatherv(*spaceComm_, latlon, latlon_global);

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

  Log::info() << "Geometry: Global tree size = " << globalTree_.size()
              << ", footprint = " << globalTree_.footprint() << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_GEOMETRY_H_
