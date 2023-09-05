/*
 * (C) Copyright 2020 UCAR.
 * (C) Crown Copyright 2023, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSERVERS_H_
#define OOPS_BASE_OBSERVERS_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Geometry.h"
#include "oops/base/GetValuePosts.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/Observer.h"
#include "oops/base/ObsOperatorBase.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/ObsTypeParameters.h"
#include "oops/base/ObsVector.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/util/Logger.h"

namespace oops {

// Note on Parameters hierarchy:
// 1. the ObserversParameters is typically the top-level Parameter for obs-related options,
//    and so is typically accessed via the "observations" key in the YAML files.
// 2. the ObserversParameters constructs an Observers object. It contains an "obsevers" key
//    the options to construct a vector of Observer objects.
template <typename MODEL, typename OBS>
class ObserversParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObserversParameters, Parameters)

  typedef typename VariableChange<MODEL>::Parameters_ VarChangeParameters_;

 public:
  Parameter<bool> obsPerturbations{"obs perturbations", false, this};
  Parameter<std::vector<ObsTypeParameters<OBS>>> observers{"observers", {}, this};
  Parameter<GetValuesParameters<MODEL>> getValues{"get values", {}, this};
};

// -----------------------------------------------------------------------------

/// \brief Computes observation operator (from GeoVaLs), applies bias correction
///        and runs QC filters
template <typename MODEL, typename OBS>
class Observers {
  typedef Geometry<MODEL>               Geometry_;
  typedef GetValues<MODEL, OBS> GetValues_;
  typedef GetValuePosts<MODEL, OBS>     GetValuePosts_;
  typedef GetValuesParameters<MODEL>    GetValuesParameters_;
  typedef ObsAuxControls<OBS>           ObsAuxCtrls_;
  typedef ObsErrors<OBS>                ObsErrors_;
  typedef Observations<OBS>             Observations_;
  typedef Observer<MODEL, OBS>          Observer_;
  typedef ObserverParameters<OBS>       ObserverParameters_;
  typedef ObsOperatorBase<OBS>          ObsOperatorBase_;
  typedef ObsSpaces<OBS>                ObsSpaces_;
  typedef ObsVector<OBS>                ObsVector_;
  typedef State<MODEL>                  State_;
  typedef PostProcessor<State_>         PostProc_;
  template <typename DATA> using ObsData_ = ObsDataVector<OBS, DATA>;
  template <typename DATA> using ObsDataVec_ = std::vector<std::shared_ptr<ObsData_<DATA>>>;

 public:
/// \brief Initializes ObsOperators, Locations, and QC data
  Observers(const ObsSpaces_ &, const std::vector<ObserverParameters_> &,
            const GetValuesParameters_ &, std::vector<std::unique_ptr<ObsOperatorBase_>>
            obsOpBases = {});
  Observers(const ObsSpaces_ &, const eckit::Configuration &,
            std::vector<std::unique_ptr<ObsOperatorBase_>> obsOpBases = {});

/// \brief Initializes variables, obs bias, obs filters (could be different for
/// different iterations
  void initialize(const Geometry_ &, const ObsAuxCtrls_ &, ObsErrors_ &,
                  PostProc_ &, const eckit::Configuration & = eckit::LocalConfiguration());

/// \brief Computes H(x) from the filled in GeoVaLs
  void finalize(Observations_ &);

  void resetObsOp(std::vector<std::unique_ptr<ObsOperatorBase_>>);

 private:
  static std::vector<ObserverParameters_> convertToParameters(const eckit::Configuration &config);
  static GetValuesParameters_ extractGetValuesParameters(const eckit::Configuration &config);

 private:
  std::vector<std::unique_ptr<Observer_>>  observers_;
  GetValuesParameters_ getValuesParams_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
Observers<MODEL, OBS>::Observers(const ObsSpaces_ & obspaces,
                                 const std::vector<ObserverParameters_> & params,
                                 const GetValuesParameters_ & getValuesParams,
                                 std::vector<std::unique_ptr<ObsOperatorBase_>> obsOpBases)
  : observers_(), getValuesParams_(getValuesParams)
{
  Log::trace() << "Observers<MODEL, OBS>::Observers start" << std::endl;
  if (obsOpBases.size() != obspaces.size()) obsOpBases.resize(obspaces.size());
  ASSERT(obspaces.size() == params.size());
  for (size_t jj = 0; jj < obspaces.size(); ++jj) {
      observers_.emplace_back(new Observer_(obspaces[jj], params[jj], std::move(obsOpBases[jj])));
  }

  Log::trace() << "Observers<MODEL, OBS>::Observers done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
Observers<MODEL, OBS>::Observers(const ObsSpaces_ & obspaces, const eckit::Configuration & config,
                                 std::vector<std::unique_ptr<ObsOperatorBase_>> obsOpBases)
  : Observers(obspaces,
              convertToParameters(config.getSubConfiguration("observers")),
              extractGetValuesParameters(config.getSubConfiguration("get values")),
              std::move(obsOpBases))
{}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void Observers<MODEL, OBS>::initialize(const Geometry_ & geom, const ObsAuxCtrls_ & obsaux,
                                       ObsErrors_ & Rmat, PostProc_ & pp,
                                       const eckit::Configuration & conf) {
  Log::trace() << "Observers<MODEL, OBS>::initialize start" << std::endl;

  std::shared_ptr<GetValuePosts_> posts(new GetValuePosts_(getValuesParams_));
  for (size_t jj = 0; jj < observers_.size(); ++jj) {
    std::vector<std::shared_ptr<GetValues_>> getvalues =
        observers_[jj]->initialize(geom, obsaux[jj], Rmat[jj], conf);
    for (std::shared_ptr<GetValues_> &gv : getvalues)
      posts->append(std::move(gv));
  }
  pp.enrollProcessor(posts);

  Log::trace() << "Observers<MODEL, OBS>::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void Observers<MODEL, OBS>::finalize(Observations_ & yobs) {
  oops::Log::trace() << "Observers<MODEL, OBS>::finalize start" << std::endl;

  for (size_t jj = 0; jj < observers_.size(); ++jj) {
    observers_[jj]->finalize(yobs[jj]);
  }

  oops::Log::trace() << "Observers<MODEL, OBS>::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void Observers<MODEL, OBS>::resetObsOp(std::vector<std::unique_ptr<ObsOperatorBase_>> obsOpBases) {
  oops::Log::trace() << "Observers<MODEL, OBS>::resetObsOp start" << std::endl;

  for (size_t jj = 0; jj < observers_.size(); ++jj) {
    observers_[jj]->resetObsOp(std::move(obsOpBases[jj]));
  }

  oops::Log::trace() << "Observers<MODEL, OBS>::resetObsOp done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
std::vector<ObserverParameters<OBS>> Observers<MODEL, OBS>::convertToParameters(
    const eckit::Configuration &config) {
  oops::Log::trace() << "Observers<MODEL, OBS>::convertToParameters start" << std::endl;

  std::vector<eckit::LocalConfiguration> subconfigs = config.getSubConfigurations();
  std::vector<ObserverParameters<OBS>> parameters(subconfigs.size());
  for (size_t i = 0; i < subconfigs.size(); ++i) {
    const eckit::LocalConfiguration &subconfig = subconfigs[i];

    // 'subconfig' will, in general, contain options irrelevant to the observer (e.g. 'obs space').
    // So we need to extract the relevant parts into a new Configuration object, 'observerConfig',
    // before validation and deserialization. Otherwise validation might fail.

    eckit::LocalConfiguration observerConfig;

    // Required keys
    observerConfig.set("obs operator", eckit::LocalConfiguration(subconfig, "obs operator"));

    // Optional keys
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

  oops::Log::trace() << "Observers<MODEL, OBS>::convertToParameters done" << std::endl;

  return parameters;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
GetValuesParameters<MODEL> Observers<MODEL, OBS>::extractGetValuesParameters(
    const eckit::Configuration & config) {
  oops::Log::trace() << "Observers<MODEL, OBS>::extractGetValuesParameters start" << std::endl;
  GetValuesParameters<MODEL> parameters{};
  parameters.deserialize(config);
  oops::Log::trace() << "Observers<MODEL, OBS>::extractGetValuesParameters done" << std::endl;
  return parameters;
}

}  // namespace oops

#endif  // OOPS_BASE_OBSERVERS_H_
