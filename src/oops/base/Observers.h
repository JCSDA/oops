/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#ifndef OOPS_BASE_OBSERVERS_H_
#define OOPS_BASE_OBSERVERS_H_

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/GetValuesPost.h"
#include "oops/base/ObsAuxControls.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsFilters.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/Locations.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

template <typename OBS>
class ObserversParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(ObserversParameters, Parameters)

 public:
  oops::RequiredParameter<eckit::LocalConfiguration> obsOperator{"obs operator", this};
  oops::Parameter<std::vector<ObsFilterParametersWrapper<OBS>>> obsFilters{"obs filters", {}, this};
};

// -----------------------------------------------------------------------------

/// \brief Computes observation operator (from GeoVaLs), applies bias correction
///        and runs QC filters
template <typename MODEL, typename OBS>
class Observers {
  typedef GeoVaLs<OBS>               GeoVaLs_;
  typedef GetValuesPost<MODEL, OBS>  GetValuesPost_;
  typedef Locations<OBS>             Locations_;
  typedef ObsAuxControls<OBS>        ObsAuxCtrls_;
  typedef ObsDiagnostics<OBS>        ObsDiags_;
  typedef Observations<OBS>          Observations_;
  typedef ObsFilters<OBS>            ObsFilters_;
  typedef ObsOperator<OBS>           ObsOperator_;
  typedef ObsSpaces<OBS>             ObsSpaces_;
  typedef std::vector<std::shared_ptr<ObsVector<OBS>>>  ObsVectors_;
  template <typename DATA> using ObsData_ = ObsDataVector<OBS, DATA>;
  template <typename DATA> using ObsDataVec_ = std::vector<std::shared_ptr<ObsData_<DATA>>>;

  typedef std::vector<std::unique_ptr<GeoVaLs_>>       GeoVaLsVec_;
  typedef std::vector<std::unique_ptr<Locations_>>     LocationsVec_;
  typedef std::vector<std::unique_ptr<ObsFilters_>>    ObsFiltersVec_;
  typedef std::vector<std::unique_ptr<ObsOperator_>>   ObsOperatorVec_;
  typedef std::vector<Variables>                       VariablesVec_;

 public:
/// \brief Initializes ObsOperators, Locations, and QC data
  Observers(const ObsSpaces_ &, const eckit::Configuration &);

/// \brief Initializes variables, obs bias, obs filters (could be different for
/// different iterations
  std::shared_ptr<PostBase<State<MODEL> > > initialize(const ObsAuxCtrls_ &,
                                                       ObsVectors_ &, const int iter = 0);

/// \brief Computes H(x) from the filled in GeoVaLs
  Observations_ finalize();

 private:
  eckit::LocalConfiguration obsconfig_;
  const ObsSpaces_ &   obspaces_;   // ObsSpaces used in H(x)
  ObsOperatorVec_      obsops_;     // Obs operators
  LocationsVec_        locations_;  // locations  (made non local by GetValuesPost)
  const ObsAuxCtrls_ * ybias_;      // Obs bias
  ObsFiltersVec_       filters_;    // QC filters
  std::shared_ptr<GetValuesPost_> getvals_;  // Postproc passed to the model during integration.
  ObsDataVec_<int>     qcflags_;    // QC flags
  int                  iterout_;    // Outer iteration
  bool                 initialized_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
Observers<MODEL, OBS>::Observers(const ObsSpaces_ & obspaces,
                                 const eckit::Configuration & config) :
  obsconfig_(config), obspaces_(obspaces), obsops_(), locations_(),
  ybias_(nullptr), filters_(), iterout_(0), initialized_(false)
{
  std::vector<eckit::LocalConfiguration> obsconfs = config.getSubConfigurations();
  for (size_t jj = 0; jj < obspaces_.size(); ++jj) {
    ObserversParameters<OBS> observerParams;
    observerParams.deserialize(obsconfs[jj]);
    /// Set up observation operators
    obsops_.emplace_back(new ObsOperator_(obspaces_[jj], observerParams.obsOperator));
    locations_.emplace_back(new Locations_(obsops_[jj]->locations()));
    qcflags_.emplace_back(new ObsData_<int>(obspaces_[jj], obspaces_[jj].obsvariables()));
  }
  Log::trace() << "Observers<MODEL, OBS> constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
std::shared_ptr<PostBase<State<MODEL> > >
Observers<MODEL, OBS>::initialize(const ObsAuxCtrls_ & obsaux, ObsVectors_ & obserrs,
                                  const int iter) {
  std::vector<eckit::LocalConfiguration> obsconfs = obsconfig_.getSubConfigurations();
  iterout_ = iter;
  ybias_ = &obsaux;
  filters_.clear();
  VariablesVec_ geovars(obspaces_.size());    // variables required from the model
  for (size_t jj = 0; jj < obspaces_.size(); ++jj) {
    ObserversParameters<OBS> observerParams;
    observerParams.deserialize(obsconfs[jj]);
    /// Set up QC filters and run preprocess
    filters_.emplace_back(new ObsFilters_(obspaces_[jj], observerParams.obsFilters,
                                          qcflags_[jj], *obserrs[jj], iterout_));
    filters_[jj]->preProcess();

    /// Set up variables requested from the model
    geovars[jj] += obsops_[jj]->requiredVars();
    geovars[jj] += (*ybias_)[jj].requiredVars();
    geovars[jj] += filters_[jj]->requiredVars();
  }

  std::vector<eckit::LocalConfiguration> gvconfs
    = util::vectoriseAndFilter(obsconfig_, "get values");

  getvals_.reset(new GetValuesPost_(obspaces_, locations_, geovars, gvconfs));

  initialized_ = true;
  Log::trace() << "Observers<MODEL, OBS>::initialize done" << std::endl;
  return getvals_;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
Observations<OBS> Observers<MODEL, OBS>::finalize() {
  oops::Log::trace() << "Observers<MODEL, OBS>::finalize start" << std::endl;
  ASSERT(initialized_);

  const GeoVaLsVec_ & geovals(getvals_->geovals());

  Observations<OBS> yobs(obspaces_);
  for (size_t jj = 0; jj < obspaces_.size(); ++jj) {
    /// call prior filters
    filters_[jj]->priorFilter(*geovals[jj]);
    /// compute H(x)
    oops::Variables vars;
    vars += filters_[jj]->requiredHdiagnostics();
    vars += (*ybias_)[jj].requiredHdiagnostics();
    ObsDiags_ ydiags(obspaces_[jj], *locations_[jj], vars);
    obsops_[jj]->simulateObs(*geovals[jj], yobs[jj], (*ybias_)[jj], ydiags);
    /// call posterior filters
    filters_[jj]->postFilter(yobs[jj], ydiags);
  }
  initialized_ = false;

// Save flags for diagnostics
  const std::string qcname  = "EffectiveQC" + std::to_string(iterout_);
  for (const auto & qcflag : qcflags_) {
    qcflag->save(qcname);
  }

  oops::Log::trace() << "Observers<MODEL, OBS>::finalize done" << std::endl;
  return yobs;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVERS_H_
