/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSSPACES_H_
#define OOPS_BASE_OBSSPACES_H_

#include <cstddef>
#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/geometry/Point2.h"

#include "oops/interface/ObsSpace.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Timer.h"

namespace oops {
  template <typename T>
  class Departures;

// -----------------------------------------------------------------------------

template <typename OBS>
class ObsSpaces : public util::Printable,
                  private util::ObjectCounter<ObsSpaces<OBS> > {
  typedef Departures<OBS>         Departures_;
  typedef ObsSpace<OBS>           ObsSpace_;

 public:
  static const std::string classname() {return "oops::ObsSpaces";}

  ObsSpaces(const eckit::Configuration &, const eckit::mpi::Comm &,
            const util::DateTime &, const util::DateTime &,
            const eckit::mpi::Comm & time = oops::mpi::myself());
  ObsSpaces(const ObsSpaces &, const eckit::geometry::Point2 &, const eckit::Configuration &);
  ~ObsSpaces();

/// Access
  std::size_t size() const {return spaces_.size();}
  ObsSpace_ & operator[](const std::size_t ii) {return *spaces_.at(ii);}
  const ObsSpace_ & operator[](const std::size_t ii) const {return *spaces_.at(ii);}

/// Assimilation window
  const util::DateTime & windowStart() const {return wbgn_;}
  const util::DateTime & windowEnd() const {return wend_;}

/// Other
  void printJo(const Departures_ &, const Departures_ &) const;  // To be changed

 private:
  void print(std::ostream &) const;

  std::vector<std::shared_ptr<ObsSpace_> > spaces_;
  const util::DateTime wbgn_;
  const util::DateTime wend_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsSpaces<OBS>::ObsSpaces(const eckit::Configuration & conf, const eckit::mpi::Comm & comm,
                          const util::DateTime & bgn, const util::DateTime & end,
                          const eckit::mpi::Comm & time)
  : spaces_(0), wbgn_(bgn), wend_(end)
{
  std::vector<eckit::LocalConfiguration> typeconfs = conf.getSubConfigurations();
  for (std::size_t jj = 0; jj < typeconfs.size(); ++jj) {
    eckit::LocalConfiguration obsconf(typeconfs[jj], "obs space");
    Log::debug() << "ObsSpaces::ObsSpaces : conf " << obsconf << std::endl;
    std::shared_ptr<ObsSpace_> tmp(new ObsSpace_(obsconf, comm, bgn, end, time));
    spaces_.push_back(tmp);
  }
  ASSERT(spaces_.size() >0);
}

// -----------------------------------------------------------------------------

template <typename OBS>
ObsSpaces<OBS>::ObsSpaces(const ObsSpaces<OBS> & obss, const eckit::geometry::Point2 & center,
                          const eckit::Configuration & conf)
  : spaces_(0), wbgn_(obss.wbgn_), wend_(obss.wend_)
{
  std::vector<eckit::LocalConfiguration> typeconfs = conf.getSubConfigurations();
  for (std::size_t jj = 0; jj < obss.size(); ++jj) {
    eckit::LocalConfiguration locconf(typeconfs[jj], "obs error.localization");
    std::shared_ptr<ObsSpace_> tmp(new ObsSpace_(obss[jj], center, locconf));
    spaces_.push_back(tmp);
  }
  ASSERT(spaces_.size() == obss.size());
}

// -----------------------------------------------------------------------------

template <typename OBS>
ObsSpaces<OBS>::~ObsSpaces() {}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsSpaces<OBS>::print(std::ostream & os) const {
  for (std::size_t jj = 0; jj < spaces_.size(); ++jj) {
    os << *spaces_[jj];
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsSpaces<OBS>::printJo(const Departures_ & dy, const Departures_ & grad) const {
  for (std::size_t jj = 0; jj < spaces_.size(); ++jj) {
    spaces_[jj]->printJo(dy[jj], grad[jj]);
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSSPACES_H_
