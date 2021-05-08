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
#include "oops/interface/GeometryIterator.h"
#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"
#include "oops/util/MemoryCounter.h"
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

// -----------------------------------------------------------------------------

template <typename OBS>
class ObsSpace : public util::Printable,
                 private util::ObjectCounter<ObsSpace<OBS> > {
  typedef typename OBS::ObsSpace  ObsSpace_;
  typedef GeometryIterator<OBS>   ObsIterator_;

 public:
  static const std::string classname() {return "oops::ObsSpace";}

  ObsSpace(const eckit::Configuration &, const eckit::mpi::Comm &,
           const util::DateTime &, const util::DateTime &,
           const eckit::mpi::Comm & time = oops::mpi::myself());
  ~ObsSpace();

  ObsSpace(const ObsSpace &) = delete;

/// Interfacing
  ObsSpace_ & obsspace() const {return *obsdb_;}  // const problem? YT

/// Assimilation window
  const util::DateTime & windowStart() const {return obsdb_->windowStart();}
  const util::DateTime & windowEnd() const {return obsdb_->windowEnd();}

  const Variables & obsvariables() const;

// Other
  const std::string & obsname() const {return obsdb_->obsname();}

  /// Iterator to the first observation
  ObsIterator_ begin() const;
  /// Iterator to the past-the-end observation (after last)
  ObsIterator_ end() const;

  const eckit::mpi::Comm & timeComm() const {return time_;}

 private:
  void print(std::ostream &) const;

  std::unique_ptr<ObsSpace_> obsdb_;
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
  util::MemoryCounter mem(classname());
  obsdb_.reset(new ObsSpace_(conf, comm, bgn, end, time));
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
GeometryIterator<OBS> ObsSpace<OBS>::begin() const {
  Log::trace() << "ObsSpace<OBS>::begin starting" << std::endl;
  util::Timer timer(classname(), "begin");
  Log::trace() << "ObsSpace<OBS>::begin done" << std::endl;
  return ObsIterator_(obsdb_->begin());
}

// -----------------------------------------------------------------------------

template <typename OBS>
GeometryIterator<OBS> ObsSpace<OBS>::end() const {
  Log::trace() << "ObsSpace<OBS>::end starting" << std::endl;
  util::Timer timer(classname(), "end");
  Log::trace() << "ObsSpace<OBS>::end done" << std::endl;
  return ObsIterator_(obsdb_->end());
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_OBSSPACE_H_
