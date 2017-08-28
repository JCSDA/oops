/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_OBSERVATIONSPACE_H_
#define OOPS_INTERFACE_OBSERVATIONSPACE_H_

#include <ostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "util/Logger.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
}

namespace oops {
  template <typename T>
  class Departures;

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObservationSpace : public util::Printable,
                         private util::ObjectCounter<ObservationSpace<MODEL> > {
  typedef Departures<MODEL>         Departures_;
  typedef typename MODEL::ObsSpace  ObsSpace_;

 public:
  static const std::string classname() {return "oops::ObservationSpace";}

  ObservationSpace(const eckit::Configuration &, const util::DateTime &, const util::DateTime &);
  ObservationSpace(const ObservationSpace &);
  ~ObservationSpace();

/// Interfacing
  ObsSpace_ & observationspace() const {return *obsdb_;}

/// Assimilation window
  const util::DateTime & windowStart() const {return obsdb_->windowStart();}
  const util::DateTime & windowEnd() const {return obsdb_->windowEnd();}
  const eckit::Configuration & config() const {return obsdb_->config();}

// Other
  void generateDistribution(const eckit::Configuration &);
  void printJo(const Departures_ &, const Departures_ &) const;

 private:
  ObservationSpace & operator=(const ObservationSpace &);
  void print(std::ostream &) const;
  boost::shared_ptr<ObsSpace_> obsdb_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ObservationSpace<MODEL>::ObservationSpace(const eckit::Configuration & conf,
                                          const util::DateTime & bgn,
                                          const util::DateTime & end) : obsdb_() {
  Log::trace() << "ObservationSpace<MODEL>::ObservationSpace starting" << std::endl;
  util::Timer timer(classname(), "ObservationSpace");
  obsdb_.reset(new ObsSpace_(conf, bgn, end));
  Log::trace() << "ObservationSpace<MODEL>::ObservationSpace done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObservationSpace<MODEL>::ObservationSpace(const ObservationSpace & other)
  : obsdb_(other.obsdb_)
{
  Log::trace() << "ObservationSpace<MODEL>::ObservationSpace copied" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObservationSpace<MODEL>::~ObservationSpace() {
  Log::trace() << "ObservationSpace<MODEL>::~ObservationSpace" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObservationSpace<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ObservationSpace<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *obsdb_;
  Log::trace() << "ObservationSpace<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObservationSpace<MODEL>::generateDistribution(const eckit::Configuration & conf) {
  Log::trace() << "ObservationSpace<MODEL>::generateDistribution starting" << std::endl;
  util::Timer timer(classname(), "generateDistribution");
  obsdb_->generateDistribution(conf);
  Log::trace() << "ObservationSpace<MODEL>::generateDistribution done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObservationSpace<MODEL>::printJo(const Departures_ & dy, const Departures_ & grad) const {
  Log::trace() << "ObservationSpace<MODEL>::printJo starting" << std::endl;
  util::Timer timer(classname(), "printJo");
  obsdb_->printJo(dy.depvalues().obsvector(), grad.depvalues().obsvector());
  Log::trace() << "ObservationSpace<MODEL>::printJo done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSERVATIONSPACE_H_
