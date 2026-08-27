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

#include "eckit/config/LocalConfiguration.h"
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

namespace detail {

// -----------------------------------------------------------------------------
/// \brief Apply global (default) settings from the section configuration
///        ("observations") to an individual obs space configuration
///
/// \details
/// For settings that appear directly in the sectionConf, apply these to the
/// individual obs space in the obsconf configuration.
///
/// A setting made on the individual obs space always takes precedence over the
/// section-level one. For now, we are using this to gloally select between the
/// ObsGroup and OSDF obs data container, but it could be used for other
/// settings in the future.
///
/// \param sectionConf the "observations" section of the configuration
/// \param obsconf an individual "obs space" configuration to be updated in place
inline void applyObsSpaceDefaults(const eckit::Configuration & sectionConf,
                                  eckit::LocalConfiguration & obsconf) {
  // ---------------------------------------------------------------------------------
  // TODO(srh): Once the migration to the OSDF container is complete, remove this block
  // of code and the "observations.obs data container" option.
  std::string container;
  if (sectionConf.get("obs data container", container) &&
      !obsconf.has("use data frame container")) {
    if (container == "OSDF") {
      obsconf.set("use data frame container", true);
    } else if (container == "ObsGroup") {
      obsconf.set("use data frame container", false);
    } else {
      throw eckit::BadValue("Unknown 'obs data container': " + container +
                            ", expected 'ObsGroup' or 'OSDF'", Here());
    }
  }
  // ---------------------------------------------------------------------------------
}

// -----------------------------------------------------------------------------

}  // namespace detail

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
  void updateObsSpaces(const eckit::Configuration &);

  /// Access
  std::size_t size() const {return spaces_.size();}
  ObsSpace_ & operator[](const std::size_t ii) {return *spaces_.at(ii);}
  const ObsSpace_ & operator[](const std::size_t ii) const {return *spaces_.at(ii);}
  bool has(const std::string &) const;

  /// Assimilation window
  const util::DateTime windowStart() const {return timeWindow_.start();}
  const util::DateTime windowEnd() const {return timeWindow_.end();}

 private:
  void print(std::ostream &) const;

  std::vector<std::shared_ptr<ObsSpace_> > spaces_;
  util::TimeWindow timeWindow_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsSpaces<OBS>::ObsSpaces(const eckit::Configuration & conf, const eckit::mpi::Comm & comm,
                          const util::TimeWindow & timeWindow,
                          const eckit::mpi::Comm & time)
  : spaces_(0), timeWindow_(timeWindow)
{
  Log::trace() << "ObsSpaces<MODEL, OBS>::ObsSpaces start" << std::endl;

  // "conf" is the whole "observations" section of the configuration, which is a mapping:
  //
  //   observations:
  //     obs data container: OSDF
  //     observers:
  //       - obs space: ...
  //       - obs space: ...
  //
  // Settings given directly under "observations:" apply to every obs space in
  // the section. Note that a bare sequence of observers is not accepted:
  // the whole section is required, so that the global settings can be seen here.
  if (!conf.has("observers")) {
    throw eckit::UserError("The 'observations' section must contain an 'observers' list. "
                           "ObsSpaces takes the whole 'observations' section, not the "
                           "'observers' list on its own.", Here());
  }
  const std::vector<eckit::LocalConfiguration> subconfigs =
      conf.getSubConfigurations("observers");

  spaces_.reserve(subconfigs.size());
  for (size_t jj = 0; jj < subconfigs.size(); ++jj) {
    eckit::LocalConfiguration obsconf(subconfigs[jj].getSubConfiguration("obs space"));
    detail::applyObsSpaceDefaults(conf, obsconf);
    auto tmp = std::make_shared<ObsSpace_>(obsconf, comm, timeWindow, time);
    spaces_.push_back(std::move(tmp));
  }
  if (spaces_.empty()) {
    Log::warning() << "ObsSpaces<MODEL, OBS>::ObsSpaces: no obs spaces created from "
                   << "configuration: " << conf << std::endl;
  }
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

// -----------------------------------------------------------------------------

template <typename OBS>
bool ObsSpaces<OBS>::has(const std::string & name) const {
  bool hasname = spaces_[0]->has(name);
  for (std::size_t jj = 1; jj < spaces_.size(); ++jj) {
    ASSERT(spaces_[jj]->has(name) == hasname);
  }
  return hasname;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsSpaces<OBS>::updateObsSpaces(const eckit::Configuration & cdaConfig) {
  Log::trace() << "ObsSpaces::appendObs start" << std::endl;
  if (cdaConfig.has("time window")) {
      util::TimeWindow newWindow(cdaConfig.getSubConfiguration("time window"));
      timeWindow_ = newWindow;
  }
  for (std::size_t jj = 0; jj < spaces_.size(); ++jj) {
    spaces_[jj]->updateObsSpace(cdaConfig);
  }
  Log::trace() << "ObsSpaces::appendObs done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSSPACES_H_
