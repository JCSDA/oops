/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSLOCALIZATIONS_H_
#define OOPS_BASE_OBSLOCALIZATIONS_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/base/Departures.h"
#include "oops/base/ObsLocalizationBase.h"
#include "oops/base/ObsSpaces.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Container for ObsLocalizations for all observation types that are used in DA
///        obsLocalizations will also loop through localizations specified for each obs space
///        and update localization factors in this loop.
template <typename MODEL, typename OBS>
class ObsLocalizations : public util::Printable,
                         private boost::noncopyable {
  typedef GeometryIterator<MODEL>  GeometryIterator_;
  typedef Departures<OBS>          Observations_;
  typedef ObsLocalizationBase<MODEL, OBS> ObsLocalization_;
  typedef ObsSpaces<OBS>           ObsSpaces_;

 public:
  static const std::string classname() {return "oops::ObsLocalizations";}

  ObsLocalizations(const eckit::Configuration &, const ObsSpaces_ &);

  void computeLocalization(const GeometryIterator_ & point,
                           Observations_ & obsvectors) const;

 private:
  void print(std::ostream &) const;
  std::vector< std::vector<std::unique_ptr<ObsLocalization_> >> local_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
ObsLocalizations<MODEL, OBS>::ObsLocalizations(const eckit::Configuration & config,
                                               const ObsSpaces_ & obspaces) {
  Log::trace() << "ObsLocalizations<MODEL, OBS>::create starting" << std::endl;
  std::vector<eckit::LocalConfiguration> obsconf = config.getSubConfigurations();
  //  loop over ob spaces
  for (size_t jj = 0; jj < obsconf.size(); ++jj) {
    //  retrieve a vector of obs localizations and loop over them
    std::vector<eckit::LocalConfiguration> obsLocConfigs =
                        obsconf[jj].getSubConfigurations("obs localizations");
    std::vector<std::unique_ptr<ObsLocalization_> > tmpVector;
    for (size_t oli = 0; oli < obsLocConfigs.size(); ++oli) {
      tmpVector.emplace_back(ObsLocalizationFactory<MODEL, OBS>::
                             create(obsLocConfigs[oli], obspaces[jj]));
    }
    //  move a vector of unique pointers form temp array to vector of vectors
    local_.emplace_back(std::move(tmpVector));
  }
  Log::trace() << "ObsLocalizations<MODEL, OBS>::create done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void ObsLocalizations<MODEL, OBS>::computeLocalization(const GeometryIterator_ & point,
                                                       Observations_ & locfactor) const {
  //  initialize locafactors to ones and then update them in the loop bellow
  locfactor.ones();
  for (size_t jj = 0; jj < local_.size(); ++jj) {
    for (size_t oli = 0; oli < local_[jj].size(); ++oli) {
      if (local_[jj][oli]) local_[jj][oli]->computeLocalization(point, locfactor[jj]);
    }
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void ObsLocalizations<MODEL, OBS>::print(std::ostream & os) const {
  for (size_t jj = 0; jj < local_.size(); ++jj) {
    for (size_t oli = 0; oli < local_[jj].size(); ++oli) {
      if (local_[jj][oli]) {
        os << "Obs space [" << jj << "] " << *local_[jj][oli];
      }
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSLOCALIZATIONS_H_
