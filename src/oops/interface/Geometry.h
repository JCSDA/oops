/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_GEOMETRY_H_
#define OOPS_INTERFACE_GEOMETRY_H_

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "oops/interface/GeometryIterator.h"
#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/ConfigurationParameter.h"
#include "oops/util/parameters/GenericParameters.h"
#include "oops/util/parameters/HasParameters_.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/ParametersOrConfiguration.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

namespace interface {
/// \brief Interface class for the geometry of the model/state space
///
/// \details Can contain information about model resolution, gridpoints, MPI distribution
///
/// Note: implementations of this interface can opt to extract their settings either from
/// a Configuration object or from a subclass of Parameters.
///
/// In the former case, they should provide a constructor with the following signature:
///
///    Geometry(const eckit::LocalConfiguration &, const eckit::mpi::Comm &);
///
/// In the latter case, the implementer should first define a subclass of Parameters
/// holding the settings of the geometry in question. The implementation of the Geometry interface
/// should then typedef `Parameters_` to the name of that subclass and provide a constructor with
/// the following signature:
///
///    Geometry(const Parameters_ &, const eckit::mpi::Comm &);
template <typename MODEL>
class Geometry : public util::Printable,
                 private util::ObjectCounter<Geometry<MODEL> > {
  typedef typename MODEL::Geometry              Geometry_;
  typedef GeometryIterator<MODEL>               GeometryIterator_;

 public:
  /// Set to Geometry_::Parameters_ if Geometry_ provides a type called Parameters_
  /// and to GenericParameters (a thin wrapper of an eckit::LocalConfiguration object) if not.
  typedef TParameters_IfAvailableElseFallbackType_t<Geometry_, GenericParameters> Parameters_;

  static const std::string classname() {return "oops::Geometry";}

  /// Constructors from yaml (and mpi communicator), implement one (using Parameters preferred)
  Geometry(const Parameters_ &, const eckit::mpi::Comm &);
  Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
  /// Constructor from pointer to the MODEL::Geometry (used in 1DVar filter)
  explicit Geometry(std::shared_ptr<const Geometry_>);
  /// Destructor (overridden for timer and log purposes)
  virtual ~Geometry();

  /// Iterator to the first gridpoint of Geometry (only used in LocalEnsembleDA)
  GeometryIterator_ begin() const;
  /// Iterator to the last gridpoint fo Geometry (only used in LocalEnsembleDA)
  GeometryIterator_ end()   const;
  /// Values of vertical coordinate in units specified by string (only used in GETKF)
  std::vector<double> verticalCoord(std::string &) const;

  /// Accessor to the geometry communicator
  const eckit::mpi::Comm & getComm() const {return geom_->getComm();}
#if ATLASIFIED
  atlas::FunctionSpace * atlasFunctionSpace() const {return geom_->atlasFunctionSpace();}
  atlas::FieldSet * atlasFieldSet() const {return geom_->atlasFieldSet();}
#endif

 protected:
  std::shared_ptr<const Geometry_> geom_;  /// pointer to the Geometry implementation

 private:
  void print(std::ostream &) const;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::Geometry(const eckit::Configuration & config,
                          const eckit::mpi::Comm & comm)
  : Geometry(validateAndDeserialize<Parameters_>(config), comm)
{}

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::Geometry(const Parameters_ & parameters,
                          const eckit::mpi::Comm & comm): geom_() {
  Log::trace() << "Geometry<MODEL>::Geometry starting" << std::endl;
  util::Timer timer(classname(), "Geometry");
  geom_.reset(new Geometry_(
                parametersOrConfiguration<HasParameters_<Geometry_>::value>(parameters),
                comm));
  Log::trace() << "Geometry<MODEL>::Geometry done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::Geometry(std::shared_ptr<const Geometry_> ptr)
  : geom_(ptr)
{
  Log::trace() << "Geometry<MODEL>::Geometry shared_ptr done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::~Geometry() {
  Log::trace() << "Geometry<MODEL>::~Geometry starting" << std::endl;
  util::Timer timer(classname(), "~Geometry");
  geom_.reset();
  Log::trace() << "Geometry<MODEL>::~Geometry done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
GeometryIterator<MODEL> Geometry<MODEL>::begin() const {
  Log::trace() << "Geometry<MODEL>::begin starting" << std::endl;
  util::Timer timer(classname(), "begin");
  Log::trace() << "Geometry<MODEL>::begin done" << std::endl;
  return GeometryIterator_(geom_->begin());
}

// -----------------------------------------------------------------------------

template <typename MODEL>
std::vector<double> Geometry<MODEL>::verticalCoord(std::string & str) const {
  Log::trace() << "Geometry<MODEL>::verticalCoord starting" << std::endl;
  util::Timer timer(classname(), "verticalCoord");
  Log::trace() << "Geometry<MODEL>::verticalCoord done" << std::endl;
  return geom_->verticalCoord(str);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
GeometryIterator<MODEL> Geometry<MODEL>::end() const {
  Log::trace() << "Geometry<MODEL>::end starting" << std::endl;
  util::Timer timer(classname(), "end");
  Log::trace() << "Geometry<MODEL>::end done" << std::endl;
  return GeometryIterator_(geom_->end());
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void Geometry<MODEL>::print(std::ostream & os) const {
  Log::trace() << "Geometry<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *geom_;
  Log::trace() << "Geometry<MODEL>::print done" << std::endl;
}

}  // namespace interface

// -----------------------------------------------------------------------------
/// \brief Geometry class used in oops; subclass of interface class above
///
/// \details Handles additional MPI communicator parameter in the constructors
/// (for MPI distribution in time, used in oops for 4DEnVar and weak-constraint 4DVar).
/// Adds extra methods that do not need to be implemented in the implementations:
/// - geometry() (for interfacing in oops)
/// - timeComm() (accessor to the MPI communicator in time)
template <typename MODEL>
class Geometry : public interface::Geometry<MODEL> {
  typedef typename MODEL::Geometry              Geometry_;
 public:
  typedef typename interface::Geometry<MODEL>::Parameters_ Parameters_;

  /// Constructor from Parameters and mpi communicators: \p geometry for spatial distribution
  /// (handled by the implementation) and \p time for distribution in time (handled by oops)
  Geometry(const Parameters_ &, const eckit::mpi::Comm & geometry,
           const eckit::mpi::Comm & time);
  /// Constructor from Configuration and mpi communicators: \p geometry for spatial distribution
  /// (handled by the implementation) and \p time for distribution in time (handled by oops)
  Geometry(const eckit::Configuration &, const eckit::mpi::Comm & geometry,
           const eckit::mpi::Comm & time = oops::mpi::myself());
  /// Constructor from pointer to the MODEL::Geometry (used in 1DVar filter)
  explicit Geometry(std::shared_ptr<const Geometry_>);

  /// Interfacing with other oops classes
  const Geometry_ & geometry() const {return *this->geom_;}

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

#endif  // OOPS_INTERFACE_GEOMETRY_H_
