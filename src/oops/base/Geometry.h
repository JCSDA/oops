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
#include "atlas/util/KDTree.h"

#include "oops/base/GeometryData.h"
#include "oops/interface/Geometry.h"
#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"
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
       cpp17::void_t<decltype(std::declval<Geometry>().closestTask(double(), double()))>>
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

  const GeometryData & generic() const {return gdata_;}

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
    return gdata_.closestTask(lat, lon);
  }
  ///@}

 private:
  const eckit::mpi::Comm * timeComm_;   /// pointer to the MPI communicator in time
  GeometryData gdata_;

  void setTrees();
};

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::Geometry(const eckit::Configuration & config,
                          const eckit::mpi::Comm & geometry, const eckit::mpi::Comm & time):
  interface::Geometry<MODEL>(config, geometry),
  timeComm_(&time),
  gdata_(this->geom_->functionSpace(), this->geom_->fields(),
         this->geom_->levelsAreTopDown(), geometry)
{
  this->setTrees();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::Geometry(const Parameters_ & parameters,
                          const eckit::mpi::Comm & geometry, const eckit::mpi::Comm & time):
  interface::Geometry<MODEL>(parameters, geometry),
  timeComm_(&time),
  gdata_(this->geom_->functionSpace(), this->geom_->fields(),
         this->geom_->levelsAreTopDown(), geometry)
{
  this->setTrees();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::Geometry(std::shared_ptr<const Geometry_> ptr):
  interface::Geometry<MODEL>(ptr),
  timeComm_(&oops::mpi::myself()),
  gdata_(this->geom_->functionSpace(), this->geom_->fields(),
         this->geom_->levelsAreTopDown(), oops::mpi::world())
{
  this->setTrees();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void Geometry<MODEL>::setTrees() {
  std::vector<double> lats;
  std::vector<double> lons;
  this->latlon(lats, lons, false);
  gdata_.setGlobalTree(lats, lons);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
bool Geometry<MODEL>::operator==(const Geometry & rhs) const {
  Log::trace() << "Geometry<MODEL>::operator== starting" << std::endl;
  bool eq = (this->geom_ == rhs.geom_ &&
             &gdata_.comm() == &rhs.gdata_.comm() && timeComm_ == rhs.timeComm_);
  Log::trace() << "Geometry<MODEL>::operator== done" << std::endl;
  return eq;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_GEOMETRY_H_
