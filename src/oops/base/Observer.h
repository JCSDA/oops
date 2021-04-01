/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#ifndef OOPS_BASE_OBSERVER_H_
#define OOPS_BASE_OBSERVER_H_

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/GetValuePost.h"
#include "oops/base/ObsFilters.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/Locations.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

template <typename OBS>
class ObserverParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(ObserverParameters, Parameters)

 public:
  oops::RequiredParameter<eckit::LocalConfiguration> obsOperator{"obs operator", this};
  oops::Parameter<std::vector<ObsFilterParametersWrapper<OBS>>> obsFilters{"obs filters", {}, this};
};

// -----------------------------------------------------------------------------

/// \brief Computes observation operator, applying bias correction and QC filters
template <typename MODEL, typename OBS>
class Observer {
  typedef Geometry<MODEL>              Geometry_;
  typedef GeoVaLs<OBS>                 GeoVaLs_;
  typedef GetValuePost<MODEL, OBS>     GetValPost_;
  typedef Locations<OBS>               Locations_;
  typedef ObsAuxControl<OBS>           ObsAuxCtrl_;
  typedef ObsDataVector<OBS, int>      ObsDataInt_;
  typedef ObsDiagnostics<OBS>          ObsDiags_;
  typedef ObsFilters<OBS>              ObsFilters_;
  typedef ObsOperator<OBS>             ObsOperator_;
  typedef ObsSpace<OBS>                ObsSpace_;
  typedef ObsVector<OBS>               ObsVector_;

 public:
/// \brief Initializes ObsOperators, Locations, and QC data
  Observer(const ObsSpace_ &, const eckit::Configuration &);

/// \brief Initializes variables, obs bias, obs filters (could be different for
/// different iterations
  std::shared_ptr<GetValPost_> initialize(const Geometry_ &, const ObsAuxCtrl_ &,
                                          ObsVector_ &, const int iter = 0);

/// \brief Computes H(x) from the filled in GeoVaLs
  void finalize(ObsVector_ &);

 private:
  eckit::LocalConfiguration     obsconfig_;
  const ObsSpace_ &             obspace_;    // ObsSpace used in H(x)
  std::unique_ptr<ObsOperator_> obsop_;      // Obs operator
  std::unique_ptr<Locations_>   locations_;  // locations
  const ObsAuxCtrl_ *           ybias_;      // Obs bias
  std::unique_ptr<ObsFilters_>  filters_;    // QC filters
  std::shared_ptr<GetValPost_>  getvals_;    // Postproc passed to the model during integration.
  std::shared_ptr<ObsDataInt_>  qcflags_;    // QC flags (should not be a pointer)
  int                           iterout_;    // Outer iteration
  bool                          initialized_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
Observer<MODEL, OBS>::Observer(const ObsSpace_ & obspace, const eckit::Configuration & config)
  : obsconfig_(config), obspace_(obspace), obsop_(), locations_(),
    ybias_(nullptr), filters_(), qcflags_(), iterout_(0), initialized_(false)
{
  Log::trace() << "Observer::Observer start" << std::endl;
  ObserverParameters<OBS> observerParams;
  observerParams.deserialize(config);
  /// Set up observation operators
  obsop_.reset(new ObsOperator_(obspace_, observerParams.obsOperator));
  qcflags_.reset(new ObsDataInt_(obspace_, obspace_.obsvariables()));

  Log::trace() << "Observer::Observer done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
std::shared_ptr<GetValuePost<MODEL, OBS>>
Observer<MODEL, OBS>::initialize(const Geometry_ & geom, const ObsAuxCtrl_ & ybias,
                                 ObsVector_ & obserr, const int iter) {
// could pass state (or even geom) and obsbias instead of control var (easier for HofX?)
  iterout_ = iter;
  ybias_ = &ybias;

  ObserverParameters<OBS> observerParams;
  observerParams.deserialize(obsconfig_);
  /// Set up QC filters and run preprocess
  filters_.reset(new ObsFilters_(obspace_, observerParams.obsFilters,
                                 qcflags_, obserr, iterout_));
  filters_->preProcess();

  locations_.reset(new Locations_(obsop_->locations()));

  /// Set up variables that will be requested from the model
  Variables geovars;
  geovars += obsop_->requiredVars();
  geovars += ybias_->requiredVars();
  geovars += filters_->requiredVars();

  eckit::LocalConfiguration gvconf = obsconfig_.getSubConfiguration("get values");

  getvals_.reset(new GetValPost_(gvconf, geom, obspace_.windowStart(),
                                 obspace_.windowEnd(), *locations_, geovars));

  initialized_ = true;
  Log::trace() << "Observer<MODEL, OBS>::initialize done" << std::endl;
  return getvals_;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void Observer<MODEL, OBS>::finalize(ObsVector_ & yobsim) {
  oops::Log::trace() << "Observer<MODEL, OBS>::finalize start" << std::endl;
  ASSERT(initialized_);

  // GetValues releases GeoVaLs, Observer takes ownership
  std::unique_ptr<GeoVaLs_> geovals = getvals_->releaseGeoVaLs();

  /// Call prior filters
  filters_->priorFilter(*geovals);

  /// Setup diagnostics
  Variables vars;
  vars += filters_->requiredHdiagnostics();
  vars += ybias_->requiredHdiagnostics();
  ObsDiags_ ydiags(obspace_, *locations_, vars);

  /// Compute H(x)
  obsop_->simulateObs(*geovals, yobsim, *ybias_, ydiags);

  /// Call posterior filters
  filters_->postFilter(yobsim, ydiags);

  /// Save flags (for diagnostics use)
  const std::string qcname  = "EffectiveQC" + std::to_string(iterout_);
  qcflags_->save(qcname);

  initialized_ = false;
  oops::Log::trace() << "Observer<MODEL, OBS>::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVER_H_
