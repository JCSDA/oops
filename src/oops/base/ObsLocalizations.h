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
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/base/Departures.h"
#include "oops/base/ObsLocalizationBase.h"
#include "oops/base/ObsSpaces.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

// -----------------------------------------------------------------------------
/// \brief Container for ObsLocalizations for all observation types that are used in DA
template <typename MODEL, typename OBS>
class ObsLocalizations : public util::Printable,
                         private boost::noncopyable {
  typedef GeometryIterator<MODEL>  GeometryIterator_;
  typedef Departures<OBS>          Observations_;
  typedef ObsLocalizationBase<MODEL, OBS> ObsLocalization_;
  typedef ObsSpaces<OBS>           ObsSpaces_;
  typedef std::vector<std::shared_ptr<ObsDataVector<OBS, int>>> ObsDataVectors_;

 public:
  static const std::string classname() {return "oops::ObsLocalizations";}

  ObsLocalizations(const eckit::Configuration &, const ObsSpaces_ &);

  /// schur-multiply \p obsvectors with observation-space localizations between
  /// observations in \p obsspaces and \p point in model-space
  void computeLocalization(const GeometryIterator_ & point,
                           ObsDataVectors_ & local, Observations_ & obsvectors) const;

 private:
  void print(std::ostream &) const;
  std::vector<std::unique_ptr<ObsLocalization_> > local_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
ObsLocalizations<MODEL, OBS>::ObsLocalizations(const eckit::Configuration & config,
                                               const ObsSpaces_ & obspaces) {
  std::vector<eckit::LocalConfiguration> obsconf = config.getSubConfigurations();
  for (size_t jj = 0; jj < obsconf.size(); ++jj) {
    eckit::LocalConfiguration conf(obsconf[jj], "obs localization");
    local_.emplace_back(ObsLocalizationFactory<MODEL, OBS>::create(conf, obspaces[jj]));
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void ObsLocalizations<MODEL, OBS>::computeLocalization(const GeometryIterator_ & point,
                            ObsDataVectors_ & local, Observations_ & obsvec) const {
  for (size_t jj = 0; jj < obsvec.size(); ++jj) {
    if (local_[jj]) local_[jj]->computeLocalization(point, *local[jj], obsvec[jj]);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void ObsLocalizations<MODEL, OBS>::print(std::ostream & os) const {
  for (size_t jj = 0; jj < local_.size(); ++jj) {
    if (local_[jj]) os << *local_[jj];
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSLOCALIZATIONS_H_
