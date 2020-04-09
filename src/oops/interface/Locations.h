/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_LOCATIONS_H_
#define OOPS_INTERFACE_LOCATIONS_H_

#include <memory>
#include <string>

#include "eckit/mpi/Comm.h"

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
class Locations : public util::Printable,
                  private util::ObjectCounter<Locations<MODEL> > {
  typedef typename MODEL::Locations             Locations_;

 public:
  static const std::string classname() {return "oops::Locations";}

  explicit Locations(const Locations_ *);
  Locations(const eckit::Configuration &, const eckit::mpi::Comm &);
  Locations(const Locations &);
  ~Locations();

/// Interfacing
  const Locations_ & locations() const {return *locs_;}

 private:
  void print(std::ostream &) const;
  std::shared_ptr<const Locations_> locs_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
Locations<MODEL>::Locations(const Locations_ * locs) : locs_(locs)
{
  Log::trace() << "Locations<MODEL>::Locations constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Locations<MODEL>::Locations(const eckit::Configuration & conf, const eckit::mpi::Comm & comm) {
  Log::trace() << "Locations<MODEL>::Locations starting" << std::endl;
  util::Timer timer(classname(), "Locations");
  locs_.reset(new Locations_(conf, comm));
  Log::trace() << "Locations<MODEL>::Locations done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Locations<MODEL>::Locations(const Locations & other)
  : locs_(other.locs_)
{
  Log::trace() << "Locations<MODEL>::Locations copied" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
Locations<MODEL>::~Locations() {
  Log::trace() << "Locations<MODEL>::~Locations starting" << std::endl;
  util::Timer timer(classname(), "~Locations");
  locs_.reset();
  Log::trace() << "Locations<MODEL>::~Locations done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Locations<MODEL>::print(std::ostream & os) const {
  Log::trace() << "Locations<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *locs_;
  Log::trace() << "Locations<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_LOCATIONS_H_
