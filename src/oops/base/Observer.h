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

#include "oops/base/Geometry.h"
#include "oops/base/GetValuePost.h"
#include "oops/base/ObsError.h"
#include "oops/base/ObsFilters.h"
#include "oops/base/Variables.h"
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

  typedef typename OBS::ObsOperator::Parameters_ ObsOperatorParameters_;

 public:
  oops::RequiredParameter<ObsOperatorParameters_> obsOperator{"obs operator", this};
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
  typedef ObsError<OBS>                ObsError_;
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
                                          ObsError_ &, const eckit::Configuration &);

/// \brief Computes H(x) from the filled in GeoVaLs
  void finalize(ObsVector_ &);

 private:
  Parameters_                   parameters_;
  const ObsSpace_ &             obspace_;    // ObsSpace used in H(x)
  std::unique_ptr<ObsOperator_> obsop_;      // Obs operator
  std::unique_ptr<Locations_>   locations_;  // locations
  const ObsAuxCtrl_ *           biascoeff_;  // bias coefficients
  ObsError_ *                   Rmat_;       // Obs error covariance
  std::unique_ptr<ObsFilters_>  filters_;    // QC filters
  std::unique_ptr<ObsVector_>   obserr_;     // Obs error std dev
  std::shared_ptr<GetValPost_>  getvals_;    // Postproc passed to the model during integration.
  std::shared_ptr<ObsDataInt_>  qcflags_;    // QC flags (should not be a pointer)
  bool                          initialized_;
  std::unique_ptr<eckit::LocalConfiguration> iterconf_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
Observer<MODEL, OBS>::Observer(const ObsSpace_ & obspace, const Parameters_ & params)
  : parameters_(params), obspace_(obspace), obsop_(), locations_(),
    biascoeff_(nullptr), filters_(), qcflags_(), initialized_(false)
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
Observer<MODEL, OBS>::initialize(const Geometry_ & geom, const ObsAuxCtrl_ & biascoeff,
                                 ObsError_ & R, const eckit::Configuration & conf) {
  Log::trace() << "Observer<MODEL, OBS>::initialize start" << std::endl;
// Save information for finalize
  iterconf_.reset(new eckit::LocalConfiguration(conf));
  biascoeff_ = &biascoeff;
  Rmat_ = &R;
  obserr_.reset(new ObsVector_(Rmat_->obserrors()));

// Set up QC filters and run preprocess
  const int iterfilt = iterconf_->getInt("iteration", 0);
  filters_.reset(new ObsFilters_(obspace_, parameters_.obsFilters,
                                 qcflags_, *obserr_, iterfilt));
  filters_->preProcess();

  locations_.reset(new Locations_(obsop_->locations()));

// Set up variables that will be requested from the model
  Variables geovars;
  geovars += obsop_->requiredVars();
  geovars += biascoeff_->requiredVars();
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
  vars += biascoeff_->requiredHdiagnostics();
  ObsDiags_ ydiags(obspace_, *locations_, vars);

  // Setup bias vector
  ObsVector_ ybias(obspace_);
  ybias.zero();

  /// Compute H(x)
  obsop_->simulateObs(*geovals, yobsim, *biascoeff_, ybias, ydiags);

  /// Call posterior filters
  filters_->postFilter(yobsim, ybias, ydiags);

  // Update R with obs errors that filters might have updated
  Rmat_->update(*obserr_);

  // Save current obs, obs error estimates and QC flags (for diagnostics use only)
  std::string siter = "";
  if (iterconf_->has("iteration")) siter = iterconf_->getString("iteration");

  if (iterconf_->getBool("save qc", true)) {
    const std::string qcname  = "EffectiveQC" + siter;
    qcflags_->save(qcname);
  }
  if (iterconf_->getBool("save hofx", true)) {
    const std::string obsname = "hofx" + siter;
    yobsim.save(obsname);
  }
  if (iterconf_->getBool("save obs errors", true)) {
    const std::string errname = "EffectiveError" + siter;
    Rmat_->save(errname);
  }
  if (iterconf_->getBool("save obs bias", true)) {
    const std::string biasname  = "ObsBias" + siter;
    ybias.save(biasname);
  }

  Log::info() << "Observer::finalize QC = " << *qcflags_ << std::endl;

  initialized_ = false;
  Log::trace() << "Observer<MODEL, OBS>::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVER_H_
