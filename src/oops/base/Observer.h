/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#ifndef OOPS_BASE_OBSERVER_H_
#define OOPS_BASE_OBSERVER_H_

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/GetValuePost.h"
#include "oops/base/ObsErrorBase.h"
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
  oops::Parameter<eckit::LocalConfiguration> getValues{
    "get values", eckit::LocalConfiguration(), this};
  oops::Parameter<eckit::LocalConfiguration> linearGetValues{
    "linear get values", eckit::LocalConfiguration(), this};
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
  typedef ObsErrorBase<OBS>            ObsError_;
  typedef ObserverParameters<OBS>      Parameters_;
  typedef ObsFilters<OBS>              ObsFilters_;
  typedef ObsOperator<OBS>             ObsOperator_;
  typedef ObsSpace<OBS>                ObsSpace_;
  typedef ObsVector<OBS>               ObsVector_;

 public:
/// \brief Initializes ObsOperators, Locations, and QC data
  Observer(const ObsSpace_ &, const Parameters_ &);

/// \brief Initializes variables, obs bias, obs filters (could be different for
/// different iterations
  std::shared_ptr<GetValPost_> initialize(const Geometry_ &, const ObsAuxCtrl_ &,
                                          ObsError_ &, const int iter);

/// \brief Computes H(x) from the filled in GeoVaLs
  void finalize(ObsVector_ &);
  void finalize(ObsVector_ &, std::shared_ptr<ObsDataInt_> &);

 private:
  Parameters_                   parameters_;
  const ObsSpace_ &             obspace_;    // ObsSpace used in H(x)
  std::unique_ptr<ObsOperator_> obsop_;      // Obs operator
  std::unique_ptr<Locations_>   locations_;  // locations
  const ObsAuxCtrl_ *           ybias_;      // Obs bias
  ObsError_ *                   Rmat_;       // Obs error covariance
  std::unique_ptr<ObsFilters_>  filters_;    // QC filters
  std::unique_ptr<ObsVector_>   obserr_;     // Obs error std dev
  std::shared_ptr<GetValPost_>  getvals_;    // Postproc passed to the model during integration.
  std::shared_ptr<ObsDataInt_>  qcflags_;    // QC flags (should not be a pointer)
  int                           iterout_;    // Outer iteration
  bool                          initialized_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
Observer<MODEL, OBS>::Observer(const ObsSpace_ & obspace, const Parameters_ & params)
  : parameters_(params), obspace_(obspace), obsop_(), locations_(),
    ybias_(nullptr), filters_(), qcflags_(), iterout_(-1), initialized_(false)
{
  Log::trace() << "Observer::Observer start" << std::endl;
  /// Set up observation operators
  obsop_.reset(new ObsOperator_(obspace_, parameters_.obsOperator));
  qcflags_.reset(new ObsDataInt_(obspace_, obspace_.obsvariables()));

  Log::trace() << "Observer::Observer done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
std::shared_ptr<GetValuePost<MODEL, OBS>>
Observer<MODEL, OBS>::initialize(const Geometry_ & geom, const ObsAuxCtrl_ & ybias,
                                 ObsError_ & R, const int iter) {
  Log::trace() << "Observer<MODEL, OBS>::initialize start" << std::endl;
// Save information for finalize
  iterout_ = iter;
  ybias_ = &ybias;
  Rmat_ = &R;
  obserr_.reset(new ObsVector_(Rmat_->obserrors()));

// Set up QC filters and run preprocess
  int iterfilt = std::max(iter, 0);
  filters_.reset(new ObsFilters_(obspace_, parameters_.obsFilters,
                                 qcflags_, *obserr_, iterfilt));
  filters_->preProcess();

  locations_.reset(new Locations_(obsop_->locations()));

// Set up variables that will be requested from the model
  Variables geovars;
  geovars += obsop_->requiredVars();
  geovars += ybias_->requiredVars();
  geovars += filters_->requiredVars();

// Set up GetValues
  getvals_.reset(new GetValPost_(parameters_.getValues, geom, obspace_.windowStart(),
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

  // Update R with obs errors that filters might have updated
  Rmat_->update(*obserr_);

  // Save current obs, obs error estimates and QC flags (for diagnostics use only)
  std::string siter = "";
  if (iterout_ >= 0) siter = std::to_string(iterout_);
  const std::string qcname  = "EffectiveQC" + siter;
  qcflags_->save(qcname);
  const std::string obsname = "hofx" + siter;
  yobsim.save(obsname);
  const std::string errname = "EffectiveError" + siter;
  Rmat_->save(errname);

  Log::info() << "Observer::finalize QC = " << *qcflags_ << std::endl;

  initialized_ = false;
  Log::trace() << "Observer<MODEL, OBS>::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void Observer<MODEL, OBS>::finalize(ObsVector_ & yobsim, std::shared_ptr<ObsDataInt_> & qc) {
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

  // Update R with obs errors that filters might have updated
  Rmat_->update(*obserr_);

  qc = qcflags_;

  Log::info() << "Observer::finalize QC = " << *qcflags_ << std::endl;

  initialized_ = false;
  Log::trace() << "Observer<MODEL, OBS>::finalize done" << std::endl;
}


// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVER_H_
