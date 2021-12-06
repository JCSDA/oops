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

#include "oops/interface/Geometry.h"
#include "oops/mpi/mpi.h"

namespace eckit {
  class Configuration;
}

namespace oops {

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

  /// Accessor to the MPI communicator for distribution in time
  const eckit::mpi::Comm & timeComm() const {return *timeComm_;}

 private:
  const eckit::mpi::Comm * timeComm_;  /// pointer to the MPI communicator in time
};

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::Geometry(const eckit::Configuration & config,
                          const eckit::mpi::Comm & geometry, const eckit::mpi::Comm & time):
  interface::Geometry<MODEL>(config, geometry), timeComm_(&time)
{}

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::Geometry(const Parameters_ & parameters,
                          const eckit::mpi::Comm & geometry, const eckit::mpi::Comm & time):
  interface::Geometry<MODEL>(parameters, geometry), timeComm_(&time)
{}

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::Geometry(std::shared_ptr<const Geometry_> ptr):
  interface::Geometry<MODEL>(ptr), timeComm_(&oops::mpi::myself())
{}

}  // namespace oops

#endif  // OOPS_BASE_GEOMETRY_H_
