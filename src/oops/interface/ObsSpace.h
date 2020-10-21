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

#include <memory>
#include <ostream>
#include <string>

#include "eckit/geometry/Point2.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
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

template <typename OBS>
class ObsSpace : public util::Printable,
                 private util::ObjectCounter<ObsSpace<OBS> > {
  typedef typename OBS::ObsSpace  ObsSpace_;
  typedef ObsVector<OBS>          ObsVector_;

 public:
  static const std::string classname() {return "oops::ObsSpace";}

  ObsSpace(const eckit::Configuration &, const eckit::mpi::Comm &,
           const util::DateTime &, const util::DateTime &,
           const eckit::mpi::Comm & time = oops::mpi::myself());
  ObsSpace(const ObsSpace &, const eckit::geometry::Point2 &,
           const eckit::Configuration &);
/// Constructor added for generic 1d-var under development in ufo
  ObsSpace(const ObsSpace_ &, const eckit::geometry::Point2 &,
           const eckit::Configuration &);
  explicit ObsSpace(const ObsSpace_ &);
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

  const eckit::mpi::Comm & timeComm() const {return time_;}

 private:
  void print(std::ostream &) const;

  std::shared_ptr<ObsSpace_> obsdb_;
  const eckit::mpi::Comm & time_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsSpace<OBS>::ObsSpace(const eckit::Configuration & conf,
                        const eckit::mpi::Comm & comm,
                        const util::DateTime & bgn,
                        const util::DateTime & end,
                        const eckit::mpi::Comm & time) : obsdb_(), time_(time) {
  Log::trace() << "ObsSpace<OBS>::ObsSpace starting" << std::endl;
  util::Timer timer(classname(), "ObsSpace");
  obsdb_.reset(new ObsSpace_(conf, comm, bgn, end, time));
  Log::trace() << "ObsSpace<OBS>::ObsSpace done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
ObsSpace<OBS>::ObsSpace(const ObsSpace<OBS> & os,
                        const eckit::geometry::Point2 & center,
                        const eckit::Configuration & conf) : obsdb_(), time_(oops::mpi::myself()) {
  Log::trace() << "ObsSpace<OBS>::ObsSpace (local) starting" << std::endl;
  util::Timer timer(classname(), "ObsSpace");
  obsdb_.reset(new ObsSpace_(os.obsspace(), center, conf));
  Log::trace() << "ObsSpace<OBS>::ObsSpace (local) done" << std::endl;
}

// -----------------------------------------------------------------------------
/// Constructor added for generic 1d-var under development in ufo
template <typename OBS>
ObsSpace<OBS>::ObsSpace(const ObsSpace_ & os, const eckit::geometry::Point2 & center,
                        const eckit::Configuration & conf): obsdb_(), time_(oops::mpi::myself()) {
  Log::trace() << "ObsSpace<OBS>::ObsSpace (local) derived state starting" << std::endl;
  util::Timer timer(classname(), "ObsSpace");
  obsdb_.reset(new ObsSpace_(os, center, conf));
  Log::trace() << "ObsSpace<OBS>::ObsSpace (local) derived state done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
ObsSpace<OBS>::ObsSpace(const ObsSpace_ & other) : obsdb_(), time_(other.time_) {
  Log::trace() << "ObsSpace<OBS>::ObsSpace starting" << std::endl;
  util::Timer timer(classname(), "ObsSpace");
  obsdb_ = other.obsdb_;
  Log::trace() << "ObsSpace<OBS>::ObsSpace done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
ObsSpace<OBS>::~ObsSpace() {
  Log::trace() << "ObsSpace<OBS>::~ObsSpace starting" << std::endl;
  util::Timer timer(classname(), "~ObsSpace");
  obsdb_.reset();
  Log::trace() << "ObsSpace<OBS>::~ObsSpace done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsSpace<OBS>::print(std::ostream & os) const {
  Log::trace() << "ObsSpace<OBS>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *obsdb_;
  Log::trace() << "ObsSpace<OBS>::print done" << std::endl;
}

// -----------------------------------------------------------------------------
//
template <typename OBS>
const Variables & ObsSpace<OBS>::obsvariables() const {
  Log::trace() << "ObsSpace<OBS>::obsvariables starting" << std::endl;
  util::Timer timer(classname(), "obsvariables");
  return obsdb_->obsvariables();
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsSpace<OBS>::printJo(const ObsVector_ & dy, const ObsVector_ & grad) const {
  Log::trace() << "ObsSpace<OBS>::printJo starting" << std::endl;
  util::Timer timer(classname(), "printJo");
  obsdb_->printJo(dy.obsvector(), grad.obsvector());
  Log::trace() << "ObsSpace<OBS>::printJo done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSSPACE_H_
