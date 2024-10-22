/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Crown Copyright 2023, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_OBSERVERSTLAD_H_
#define OOPS_BASE_OBSERVERSTLAD_H_

#include <algorithm>
#include <memory>
#include <utility>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/Departures.h"
#include "oops/base/Geometry.h"
#include "oops/base/GetValueTLADs.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsAuxIncrements.h"
#include "oops/base/Observers.h"
#include "oops/base/ObserverTLAD.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/util/DateTime.h"

namespace oops {

/// Computes observation equivalent TL and AD to/from increments.

template <typename MODEL, typename OBS>
class ObserversTLAD {
  typedef Departures<OBS>             Departures_;
  typedef Geometry<MODEL>             Geometry_;
  typedef GeoVaLs<OBS>                GeoVaLs_;
  typedef GetValues<MODEL, OBS>       GetValues_;
  typedef GetValueTLADs<MODEL, OBS>   GetValueTLADs_;
  typedef GetValuesParameters<MODEL>  GetValuesParameters_;
  typedef Observations<OBS>           Observations_;
  typedef ObsAuxControls<OBS>         ObsAuxCtrls_;
  typedef ObsAuxIncrements<OBS>       ObsAuxIncrs_;
  typedef ObsDataVector<OBS, int>     ObsDataInt_;
  typedef ObserverParameters<OBS>     ObserverParameters_;
  typedef ObserverTLAD<MODEL, OBS>    ObserverTLAD_;
  typedef ObsSpaces<OBS>              ObsSpaces_;
  typedef PostProcessorTLAD<MODEL>    PostProcTLAD_;

 public:
  ObserversTLAD(const ObsSpaces_ &, const std::vector<ObserverParameters<OBS>> &,
                const GetValuesParameters_ &);
  ObserversTLAD(const ObsSpaces_ &, const eckit::Configuration &);

  void initializeTraj(const Geometry_ &, const ObsAuxCtrls_ &, PostProcTLAD_ &);
  void finalizeTraj(const std::vector<ObsDataInt_> &);

  void initializeTL(PostProcTLAD_ &);
  void finalizeTL(const ObsAuxIncrs_ &, Departures_ &);

  void initializeAD(const Departures_ &, ObsAuxIncrs_ &, PostProcTLAD_ &);
  void finalizeAD();

/// Accessors
  ObserverTLAD_ & operator[](const std::size_t ii) {return *observers_.at(ii);}
  const std::shared_ptr<GetValueTLADs_> & posts() {return posts_;}

 private:
  static std::vector<ObserverParameters_> convertToParameters(const eckit::Configuration &config);
  static GetValuesParameters_ extractGetValuesParameters(const eckit::Configuration & config);

  std::vector<std::unique_ptr<ObserverTLAD_>>  observers_;
  std::shared_ptr<GetValueTLADs_> posts_;
  util::DateTime winbgn_;
  util::DateTime winend_;
  GetValuesParameters_ getValuesParams_;
};

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
ObserversTLAD<MODEL, OBS>::ObserversTLAD(const ObsSpaces_ & obspaces,
                                         const std::vector<ObserverParameters<OBS>> & obsParams,
                                         const GetValuesParameters_ & getValuesParams)
  : observers_(), winbgn_(obspaces.windowStart()), winend_(obspaces.windowEnd()),
    getValuesParams_(getValuesParams)
{
  Log::trace() << "ObserversTLAD<MODEL, OBS>::ObserversTLAD start" << std::endl;
  for (size_t jj = 0; jj < obspaces.size(); ++jj) {
    const bool passive = obsParams[jj].monitoringOnly;
    std::unique_ptr<ObserverTLAD_> tmp;
    if (!passive) tmp.reset(new ObserverTLAD_(obspaces[jj], obsParams[jj]));
    observers_.push_back(std::move(tmp));
  }
  Log::trace() << "ObserversTLAD<MODEL, OBS>::ObserversTLAD done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
ObserversTLAD<MODEL, OBS>::ObserversTLAD(const ObsSpaces_ & obspaces,
                                         const eckit::Configuration & config)
  : ObserversTLAD(obspaces,
                  convertToParameters(config.getSubConfiguration("observers")),
                  extractGetValuesParameters(config.getSubConfiguration("get values")))
//  : ObserversTLAD(obspaces, convertToParameters(config.getSubConfiguration("observers")))
{}

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserversTLAD<MODEL, OBS>::initializeTraj(const Geometry_ & geom, const ObsAuxCtrls_ & ybias,
                                               PostProcTLAD_ & pp) {
  Log::trace() << "ObserversTLAD<MODEL, OBS>::initializeTraj start" << std::endl;
  posts_.reset(new GetValueTLADs_(getValuesParams_, winbgn_, winend_));
  for (size_t jj = 0; jj < observers_.size(); ++jj) {
    if (observers_[jj]) {
      for (std::shared_ptr<GetValues_> &getvalues : observers_[jj]->initializeTraj(geom, ybias[jj]))
        posts_->append(std::move(getvalues));
    }
  }
  pp.enrollProcessor(posts_);
  Log::trace() << "ObserversTLAD<MODEL, OBS>::initializeTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserversTLAD<MODEL, OBS>::finalizeTraj(const std::vector<ObsDataInt_> & qcflags) {
  Log::trace() << "ObserversTLAD<MODEL, OBS>::finalizeTraj start" << std::endl;
  for (size_t jj = 0; jj < observers_.size(); ++jj) {
    if (observers_[jj]) observers_[jj]->finalizeTraj(qcflags[jj]);
  }
  Log::trace() << "ObserversTLAD<MODEL, OBS>::finalizeTraj done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserversTLAD<MODEL, OBS>::initializeTL(PostProcTLAD_ & pp) {
  Log::trace() << "ObserversTLAD<MODEL, OBS>::initializeTL start" << std::endl;
  pp.enrollProcessor(posts_);
  Log::trace() << "ObserversTLAD<MODEL, OBS>::initializeTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserversTLAD<MODEL, OBS>::finalizeTL(const ObsAuxIncrs_ & ybias, Departures_ & dy) {
  Log::trace() << "ObserversTLAD<MODEL, OBS>::finalizeTL start" << std::endl;
  for (size_t jj = 0; jj < observers_.size(); ++jj) {
    if (observers_[jj]) observers_[jj]->finalizeTL(ybias[jj], dy[jj]);
  }
  Log::trace() << "ObserversTLAD<MODEL, OBS>::finalizeTL done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserversTLAD<MODEL, OBS>::initializeAD(const Departures_ & dy, ObsAuxIncrs_ & ybias,
                                             PostProcTLAD_ & pp) {
  Log::trace() << "ObserversTLAD<MODEL, OBS>::initializeAD start" << std::endl;
  for (size_t jj = 0; jj < observers_.size(); ++jj) {
    if (observers_[jj]) observers_[jj]->initializeAD(dy[jj], ybias[jj]);
  }
  pp.enrollProcessor(posts_);
  Log::trace() << "ObserversTLAD<MODEL, OBS>::initializeAD done" << std::endl;
}
// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
void ObserversTLAD<MODEL, OBS>::finalizeAD() {
  Log::trace() << "ObserversTLAD<MODEL, OBS>::finalizeAD start" << std::endl;
  for (size_t jj = 0; jj < observers_.size(); ++jj) {
    if (observers_[jj]) observers_[jj]->finalizeAD();
  }
  Log::trace() << "ObserversTLAD<MODEL, OBS>::finalizeAD done" << std::endl;
}
// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
std::vector<ObserverParameters<OBS>> ObserversTLAD<MODEL, OBS>::convertToParameters(
    const eckit::Configuration &config) {
  Log::trace() << "ObserversTLAD<MODEL, OBS>::convertToParameters start" << std::endl;
  Log::trace() << "ObserversTLAD::convertToParameters conf " << config << std::endl;

  std::vector<eckit::LocalConfiguration> subconfigs = config.getSubConfigurations();
  Log::trace() << "ObserversTLAD::convertToParameters size " << subconfigs.size() << std::endl;
  std::vector<ObserverParameters<OBS>> parameters(subconfigs.size());
  for (size_t i = 0; i < subconfigs.size(); ++i) {
    const eckit::LocalConfiguration &subconfig = subconfigs[i];

    // 'subconfig' will, in general, contain options irrelevant to the observer (e.g. 'obs space').
    // So we need to extract the relevant parts into a new Configuration object, 'observerConfig',
    // before validation and deserialization. Otherwise validation might fail.

    eckit::LocalConfiguration observerConfig;

    // Required keys
    observerConfig.set("obs operator", eckit::LocalConfiguration(subconfig, "obs operator"));
    observerConfig.set("monitoring only", subconfig.getBool("monitoring only", false));

    if (subconfig.has("linear obs operator"))
      observerConfig.set("linear obs operator",
                         eckit::LocalConfiguration(subconfig, "linear obs operator"));

    // Optional keys
    eckit::LocalConfiguration filterConfig;
    if (subconfig.get("obs filtering", filterConfig))
      observerConfig.set("obs filtering", filterConfig);
    std::vector<eckit::LocalConfiguration> filterConfigs;
    if (subconfig.get("obs filters", filterConfigs))
      observerConfig.set("obs filters", filterConfigs);
    if (subconfig.get("obs pre filters", filterConfigs))
      observerConfig.set("obs pre filters", filterConfigs);
    if (subconfig.get("obs prior filters", filterConfigs))
      observerConfig.set("obs prior filters", filterConfigs);
    if (subconfig.get("obs post filters", filterConfigs))
      observerConfig.set("obs post filters", filterConfigs);
    eckit::LocalConfiguration getValuesConfig;
    if (subconfig.get("get values", getValuesConfig))
      observerConfig.set("get values", getValuesConfig);

    parameters[i].deserialize(observerConfig);
  }

  Log::trace() << "ObserversTLAD<MODEL, OBS>::convertToParameters done" << std::endl;

  return parameters;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
GetValuesParameters<MODEL> ObserversTLAD<MODEL, OBS>::extractGetValuesParameters(
    const eckit::Configuration & config) {
  GetValuesParameters<MODEL> parameters{};
  parameters.deserialize(config);
  return parameters;
}

}  // namespace oops

#endif  // OOPS_BASE_OBSERVERSTLAD_H_
