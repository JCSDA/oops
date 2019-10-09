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

#include "eckit/mpi/Comm.h"

#include "oops/interface/GeometryIterator.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class Geometry : public util::Printable,
                 private util::ObjectCounter<Geometry<MODEL> > {
  typedef typename MODEL::Geometry              Geometry_;
  typedef GeometryIterator<MODEL>               GeometryIterator_;

 public:
  static const std::string classname() {return "oops::Geometry";}

  Geometry(const eckit::Configuration &, const eckit::mpi::Comm &);
  Geometry(const Geometry &);
  explicit Geometry(boost::shared_ptr<const Geometry_>);
  ~Geometry();

/// Interfacing
  const Geometry_ & geometry() const {return *geom_;}

  GeometryIterator_ begin() const;
  GeometryIterator_ end()   const;

  const eckit::mpi::Comm & getComm() const {return geom_->getComm();}

 private:
  Geometry & operator=(const Geometry &);
  void print(std::ostream &) const;
  boost::shared_ptr<const Geometry_> geom_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::Geometry(const eckit::Configuration & conf,
                          const eckit::mpi::Comm & comm): geom_() {
  Log::trace() << "Geometry<MODEL>::Geometry starting" << std::endl;
  util::Timer timer(classname(), "Geometry");
  geom_.reset(new Geometry_(conf, comm));
  Log::trace() << "Geometry<MODEL>::Geometry done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::Geometry(const Geometry & other): geom_(other.geom_) {
  Log::trace() << "Geometry<MODEL>::Geometry copy done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Geometry<MODEL>::Geometry(boost::shared_ptr<const Geometry_> ptr): geom_(ptr) {
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
