/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_OBSSPACE_H_
#define OOPS_INTERFACE_OBSSPACE_H_

#include <ostream>
#include <string>

#include <boost/shared_ptr.hpp>

#include "eckit/geometry/Point2.h"
#include "eckit/mpi/Comm.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace util {
  class DateTime;
}

namespace oops {
  template <typename T> class ObsVector;

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsSpace : public util::Printable,
                 private boost::noncopyable,
                 private util::ObjectCounter<ObsSpace<MODEL> > {
  typedef typename MODEL::ObsSpace  ObsSpace_;
  typedef ObsVector<MODEL>          ObsVector_;

 public:
  static const std::string classname() {return "oops::ObsSpace";}

  ObsSpace(const eckit::Configuration &, const eckit::mpi::Comm &,
           const util::DateTime &, const util::DateTime &);
  ObsSpace(const ObsSpace &, const eckit::geometry::Point2 &,
           const double &, const int &);
  ObsSpace(const ObsSpace_ &, const eckit::geometry::Point2 &,
           const double &, const int &);
  ~ObsSpace();

/// Interfacing
  ObsSpace_ & obsspace() const {return *obsdb_;}  // const problem? YT

/// Assimilation window
  const util::DateTime & windowStart() const {return obsdb_->windowStart();}
  const util::DateTime & windowEnd() const {return obsdb_->windowEnd();}

  const Variables & obsvariables() const;

// Other
  void printJo(const ObsVector_ &, const ObsVector_ &) const;
  const std::string & obsname() const {return obsdb_->obsname();}

 private:
  void print(std::ostream &) const;
  boost::shared_ptr<ObsSpace_> obsdb_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsSpace<MODEL>::ObsSpace(const eckit::Configuration & conf,
                          const eckit::mpi::Comm & comm,
                          const util::DateTime & bgn,
                          const util::DateTime & end) : obsdb_() {
  Log::trace() << "ObsSpace<MODEL>::ObsSpace starting" << std::endl;
  util::Timer timer(classname(), "ObsSpace");
  obsdb_.reset(new ObsSpace_(conf, comm, bgn, end));
  Log::trace() << "ObsSpace<MODEL>::ObsSpace done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsSpace<MODEL>::ObsSpace(const ObsSpace<MODEL> & os,
                          const eckit::geometry::Point2 & center,
                          const double & dist, const int & maxnum) : obsdb_() {
  Log::trace() << "ObsSpace<MODEL>::ObsSpace (local) starting" << std::endl;
  util::Timer timer(classname(), "ObsSpace");
  obsdb_.reset(new ObsSpace_(os.obsspace(), center, dist, maxnum));
  Log::trace() << "ObsSpace<MODEL>::ObsSpace (local) done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsSpace<MODEL>::ObsSpace(const ObsSpace_ & os,
     const eckit::geometry::Point2 & center, const double & dist, const int & maxnum):
         obsdb_() {
  Log::trace() << "ObsSpace<MODEL>::ObsSpace (local) derived state starting" << std::endl;
  util::Timer timer(classname(), "ObsSpace");
  obsdb_.reset(new ObsSpace_(os, center, dist, maxnum));
  Log::trace() << "ObsSpace<MODEL>::ObsSpace (local) derived state done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsSpace<MODEL>::~ObsSpace() {
  Log::trace() << "ObsSpace<MODEL>::~ObsSpace starting" << std::endl;
  util::Timer timer(classname(), "~ObsSpace");
  obsdb_.reset();
  Log::trace() << "ObsSpace<MODEL>::~ObsSpace done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObsSpace<MODEL>::print(std::ostream & os) const {
  Log::trace() << "ObsSpace<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *obsdb_;
  Log::trace() << "ObsSpace<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------
//
template <typename MODEL>
const Variables & ObsSpace<MODEL>::obsvariables() const {
  Log::trace() << "ObsSpace<MODEL>::obsvariables starting" << std::endl;
  util::Timer timer(classname(), "obsvariables");
  return obsdb_->obsvariables();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObsSpace<MODEL>::printJo(const ObsVector_ & dy, const ObsVector_ & grad) const {
  Log::trace() << "ObsSpace<MODEL>::printJo starting" << std::endl;
  util::Timer timer(classname(), "printJo");
  obsdb_->printJo(dy.obsvector(), grad.obsvector());
  Log::trace() << "ObsSpace<MODEL>::printJo done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSSPACE_H_
