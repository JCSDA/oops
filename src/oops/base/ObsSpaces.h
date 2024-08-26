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
#include <utility>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/interface/ObsSpace.h"
#include "oops/mpi/mpi.h"
#include "oops/util/ConfigFunctions.h"  // for vectoriseAndFilter
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/TimeWindow.h"

namespace oops {

// -----------------------------------------------------------------------------
template <typename OBS>
class ObsSpaces : public util::Printable,
                  private util::ObjectCounter<ObsSpaces<OBS> > {
  typedef ObsSpace<OBS>                   ObsSpace_;

 public:
  static const std::string classname() {return "oops::ObsSpaces";}

  ObsSpaces(const eckit::Configuration &, const eckit::mpi::Comm &,
            const util::TimeWindow &,
            const eckit::mpi::Comm & time = oops::mpi::myself());
  ~ObsSpaces();

/// Save files
  void save() const;

/// Append new obs
  void appendObs(const eckit::Configuration & appendConfig);

/// Access
  std::size_t size() const {return spaces_.size();}
  ObsSpace_ & operator[](const std::size_t ii) {return *spaces_.at(ii);}
  const ObsSpace_ & operator[](const std::size_t ii) const {return *spaces_.at(ii);}

/// Assimilation window
const util::DateTime windowStart() const {return timeWindow_.start();}
const util::DateTime windowEnd() const {return timeWindow_.end();}

 private:
  void print(std::ostream &) const;

  std::vector<std::shared_ptr<ObsSpace_> > spaces_;
  const util::TimeWindow timeWindow_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsSpaces<OBS>::ObsSpaces(const eckit::Configuration & conf, const eckit::mpi::Comm & comm,
                          const util::TimeWindow & timeWindow,
                          const eckit::mpi::Comm & time)
  : spaces_(0), timeWindow_(timeWindow)
{
  Log::trace() << "ObsSpaces<MODEL, OBS>::ObsSpaces start" << std::endl;
  std::vector<eckit::LocalConfiguration> subconfigs = conf.getSubConfigurations();
  spaces_.reserve(subconfigs.size());
  for (size_t jj = 0; jj < subconfigs.size(); ++jj) {
    const eckit::LocalConfiguration obsconf(subconfigs[jj].getSubConfiguration("obs space"));
    auto tmp = std::make_shared<ObsSpace_>(obsconf, comm, timeWindow, time);
    spaces_.push_back(std::move(tmp));
  }
  ASSERT(spaces_.size() >0);
  Log::trace() << "ObsSpaces<MODEL, OBS>::ObsSpaces done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
ObsSpaces<OBS>::~ObsSpaces() {}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsSpaces<OBS>::save() const {
  for (std::size_t jj = 0; jj < spaces_.size(); ++jj) {
    spaces_[jj]->save();
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsSpaces<OBS>::print(std::ostream & os) const {
  for (std::size_t jj = 0; jj < spaces_.size(); ++jj) {
    os << *spaces_[jj];
  }
}

template<typename OBS>
void ObsSpaces<OBS>::appendObs(const eckit::Configuration & appendConfig) {
  Log::trace() << "ObsSpaces::appendObs start" << std::endl;
  std::string appendDir = appendConfig.getString("obs append directory");
  for (std::size_t jj = 0; jj < spaces_.size(); ++jj) {
    spaces_[jj]->append(appendDir);
  }
  Log::trace() << "ObsSpaces::appendObs done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSSPACES_H_
