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

#include <string>

#include <boost/scoped_ptr.hpp>

#include "util/Logger.h"
#include "oops/interface/ObservationSpace.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "util/Timer.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class Locations : public util::Printable,
                  private util::ObjectCounter<Locations<MODEL> > {
  typedef typename MODEL::Locations             Locations_;

 public:
  static const std::string classname() {return "oops::Locations";}

  Locations(const ObservationSpace<MODEL> &,
            const util::DateTime &, const util::DateTime &);
  ~Locations();

/// Interfacing
  const Locations_ & locations() const {return *locs_;}

 private:
  void print(std::ostream &) const;
  boost::scoped_ptr<Locations_> locs_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
Locations<MODEL>::Locations(const ObservationSpace<MODEL> & os,
                            const util::DateTime & t1, const util::DateTime & t2)
  : locs_()
{
  Log::trace() << "Locations<MODEL>::Locations starting" << std::endl;
  util::Timer timer(classname(), "Locations");
  locs_.reset(new Locations_(os.observationspace(), t1, t2));
  Log::trace() << "Locations<MODEL>::Locations done" << std::endl;
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
//  os << *increment_;
  Log::trace() << "Locations<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_LOCATIONS_H_
