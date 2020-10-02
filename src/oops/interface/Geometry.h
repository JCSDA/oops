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

#include <ostream>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

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

// -----------------------------------------------------------------------------

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

  Geometry(const Parameters_ &, const eckit::mpi::Comm &, const eckit::mpi::Comm &);
  Geometry(const eckit::Configuration &, const eckit::mpi::Comm &,
           const eckit::mpi::Comm & time = oops::mpi::myself());
  Geometry(const Geometry &);
  explicit Geometry(boost::shared_ptr<const Geometry_>);
  ~Geometry();

/// Interfacing
  const Geometry_ & geometry() const {return *geom_;}

  // begining and end of the grid counter on this mpi tile
  GeometryIterator_ begin() const;
  GeometryIterator_ end()   const;
  // vertical coordinate in units specified by string
  std::vector<double> verticalCoord(std::string &) const;

  const eckit::mpi::Comm & getComm() const {return geom_->getComm();}
#if ATLASIFIED
  atlas::FunctionSpace * atlasFunctionSpace() const {return geom_->atlasFunctionSpace();}
  atlas::FieldSet * atlasFieldSet() const {return geom_->atlasFieldSet();}
#endif

  const eckit::mpi::Comm & timeComm() const {return time_;}
 private:
  Geometry & operator=(const Geometry &);
  void print(std::ostream &) const;
  boost::shared_ptr<const Geometry_> geom_;
  const eckit::mpi::Comm & time_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::Geometry(const eckit::Configuration & config,
                          const eckit::mpi::Comm & comm,
                          const eckit::mpi::Comm & time)
  : Geometry(validateAndDeserialize<Parameters_>(config), comm, time)
{}

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::Geometry(const Parameters_ & parameters,
                          const eckit::mpi::Comm & comm,
                          const eckit::mpi::Comm & time): geom_(), time_(time) {
  Log::trace() << "Geometry<MODEL>::Geometry starting" << std::endl;
  util::Timer timer(classname(), "Geometry");
  geom_.reset(new Geometry_(
                parametersOrConfiguration<HasParameters_<Geometry_>::value>(parameters),
                comm));
  Log::trace() << "Geometry<MODEL>::Geometry done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::Geometry(const Geometry & other): geom_(other.geom_), time_(other.time_) {
  Log::trace() << "Geometry<MODEL>::Geometry copy done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::Geometry(boost::shared_ptr<const Geometry_> ptr)
  : geom_(ptr), time_(oops::mpi::myself())
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

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_GEOMETRY_H_
