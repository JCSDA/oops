/*
 * (C) Copyright 2020 UCAR.
 * (C) Crown Copyright 2023, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#ifndef OOPS_BASE_OBSERVER_H_
#define OOPS_BASE_OBSERVER_H_

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Geometry.h"
#include "oops/base/GetValues.h"
#include "oops/base/GetValueTLADs.h"
#include "oops/base/Locations.h"
#include "oops/base/ObsOperatorBase.h"
#include "oops/base/ObsVariables.h"
#include "oops/base/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsError.h"
#include "oops/interface/ObsFilter.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Computes observation operator, applying bias correction and QC filter
template <typename MODEL, typename OBS>
class Observer {
  typedef Geometry<MODEL>              Geometry_;
  typedef GeoVaLs<OBS>                 GeoVaLs_;
  typedef GetValueTLADs<MODEL, OBS>    GetValueTLADs_;
  typedef GetValues<MODEL, OBS>        GetValues_;
  typedef Locations<OBS>               Locations_;
  typedef ObsAuxControl<OBS>           ObsAuxCtrl_;
  typedef ObsDataVector<OBS, int>      ObsDataInt_;
  typedef ObsDiagnostics<OBS>          ObsDiags_;
  typedef ObsError<OBS>                ObsError_;
  typedef ObsFilter<OBS>               ObsFilter_;
  typedef ObsOperator<OBS>             ObsOperator_;
  typedef ObsOperatorBase<OBS>         ObsOperatorBase_;
  typedef ObsSpace<OBS>                ObsSpace_;
  typedef ObsVector<OBS>               ObsVector_;
  typedef ObsDataVector<OBS, float>    ObsDataVector_;

 public:
/// \brief Initializes ObsOperators, Locations, and QC data
  Observer(const ObsSpace_ & obspace, const eckit::Configuration & conf,
           std::unique_ptr<ObsOperatorBase_> obsOpBase = nullptr);

/// \brief Initializes variables, obs bias, obs filter (could be different for
/// different iterations
  std::shared_ptr<GetValues_> initialize(const Geometry_ &, const ObsAuxCtrl_ &,
                                         ObsError_ &, const eckit::Configuration &);

/// \brief Computes H(x) from the filled in GeoVaLs
  void finalize(ObsVector_ &, ObsDataInt_ &);

  void resetObsPert(const Geometry_ &, std::unique_ptr<ObsOperatorBase_>,
                    const std::shared_ptr<GetValues_> &);

  void updateObserver(const eckit::Configuration &);

 private:
  const ObsSpace_ &                 obspace_;       // ObsSpace used in H(x)
  Variables                         geovars_;       // All required variables
  std::vector<size_t>               varsizes_;      // Sizes of these variables
  std::unique_ptr<ObsOperatorBase_> obsop_;         // Obs operator
  std::unique_ptr<Locations_>       locations_;     // Obs locations
  const ObsAuxCtrl_ *               biascoeff_;     // bias coefficients
  ObsError_ *                       Rmat_;          // Obs error covariance
  std::unique_ptr<ObsFilter_>       filter_;        // QC filter
  std::shared_ptr<ObsDataVector_>   obserrfilter_;  // Obs error std dev for processed variables
  // Instances of GetValues. Each receives a list of model variables and a set of paths along which
  // these variables should be interpolated. The interpolated values are stored in a single GeoVaLs
  // object (shared between all instances of GetValues).
  std::shared_ptr<GetValues_>       getvals_;
  std::shared_ptr<ObsDataInt_>      qcflags_;       // QC flags (should not be a pointer)
  bool                              initialized_;
  bool                              hasGeoVaLsFile_;
  std::unique_ptr<eckit::LocalConfiguration> iterconf_;
  eckit::LocalConfiguration gvConf_;
  eckit::LocalConfiguration filterConf_;
  eckit::LocalConfiguration geovalsConf_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
Observer<MODEL, OBS>::Observer(const ObsSpace_ & obspace, const eckit::Configuration & conf,
                               std::unique_ptr<ObsOperatorBase_> obsOpBase)
  : obspace_(obspace), geovars_(), varsizes_(), obsop_(), biascoeff_(nullptr), filter_(),
    qcflags_(), initialized_(false), hasGeoVaLsFile_(false),
    gvConf_(conf.getSubConfiguration("get values")), filterConf_(), geovalsConf_()
{
  Log::trace() << "Observer::Observer start" << std::endl;

  /// Set up observation operators
  if (obsOpBase == nullptr) {
    obsop_.reset(new ObsOperator_(obspace_, eckit::LocalConfiguration(conf, "obs operator")));
  } else {
    obsop_ = std::move(obsOpBase);
  }
  qcflags_.reset(new ObsDataInt_(obspace_, obspace_.obsvariables()));
  obserrfilter_.reset(new ObsDataVector_(obspace_, obspace_.obsvariables(), "ObsError"));
  eckit::LocalConfiguration tmpconf;
  if (conf.get("obs filtering", tmpconf))     filterConf_.set("obs filtering", tmpconf);
  if (conf.get("obs filters", tmpconf))       filterConf_.set("obs filters", tmpconf);
  if (conf.get("obs pre filters", tmpconf))   filterConf_.set("obs pre filters", tmpconf);
  if (conf.get("obs prior filters", tmpconf)) filterConf_.set("obs prior filters", tmpconf);
  if (conf.get("obs post filters", tmpconf))  filterConf_.set("obs post filters", tmpconf);

  if (conf.get("geovals", geovalsConf_)) {
    hasGeoVaLsFile_ = true;
  }
  Log::trace() << "Observer::Observer done" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
std::shared_ptr<GetValues<MODEL, OBS>>
Observer<MODEL, OBS>::initialize(const Geometry_ & geom, const ObsAuxCtrl_ & biascoeff,
                                 ObsError_ & R, const eckit::Configuration & conf) {
  Log::trace() << "Observer<MODEL, OBS>::initialize start" << std::endl;
// Save information for finalize
  iterconf_.reset(new eckit::LocalConfiguration(conf));
  biascoeff_ = &biascoeff;
  Rmat_ = &R;

  // Set up QC filter and run preprocess
  const int iterfilt = iterconf_->getInt("iteration", 0);
  filter_.reset(new ObsFilter_(obspace_, filterConf_, qcflags_, obserrfilter_, iterfilt));
  filter_->preProcess();

  if (!initialized_) {
// Get the list of required variables
    geovars_ = obsop_->requiredVars();
    geovars_ += biascoeff_->requiredVars();
    geovars_ += filter_->requiredVars();
    varsizes_ = geom.variableSizes(geovars_);

// Get the observation locations and their discretizations
    locations_ = std::make_unique<Locations_>(obsop_->locations());

// Set up GetValues
    if (!hasGeoVaLsFile_) {
      getvals_.reset(new GetValues_(gvConf_, geom, obspace_.timeWindow(),
                                    *locations_, geovars_));
    }
    initialized_ = true;
  }

  Log::trace() << "Observer<MODEL, OBS>::initialize done" << std::endl;
  return getvals_;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void Observer<MODEL, OBS>::finalize(ObsVector_ & yobsim, ObsDataInt_ & qcflags) {
  oops::Log::trace() << "Observer<MODEL, OBS>::finalize start" << std::endl;
  ASSERT(initialized_);

  // Fill GeoVaLs
  GeoVaLs_ geovals = hasGeoVaLsFile_
                     ? GeoVaLs_(geovalsConf_, obspace_, geovars_)
                     : GeoVaLs_(*locations_, geovars_, varsizes_);
  if (!hasGeoVaLsFile_) {
    if (getvals_->useMethodsTL()) {
      getvals_->fillGeoVaLsTL(geovals);
    } else {
      getvals_->fillGeoVaLs(geovals);
    }
  }

  // Compute the reduced representation of the GeoVaLs for which it's been requested
  oops::Variables reducedVars = biascoeff_->requiredVars();
  reducedVars += filter_->requiredVars();
  obsop_->computeReducedVars(reducedVars, geovals);

  /// Call prior filter
  filter_->priorFilter(geovals);

  /// Setup diagnostics
  ObsVariables vars;
  vars += filter_->requiredHdiagnostics();
  vars += biascoeff_->requiredHdiagnostics();
  // The current interface makes it possible to assign different location sampling methods not only
  // to GeoVaLs, but also to ObsDiagnostics. We could simplify things and assume there'll always
  // be a 1-to-1 mapping between obs locations and columns of ObsDiagnostics.
  ObsDiags_ ydiags(obspace_, *locations_, vars);

  // Setup bias vector
  ObsVector_ ybias(obspace_);
  ybias.zero();

  /// Compute H(x)
  obsop_->simulateObs(geovals, yobsim, *biascoeff_, *qcflags_, ybias, ydiags);

  /// Call posterior filter
  filter_->postFilter(geovals, yobsim, ybias, ydiags);
  obserrfilter_->mask(*qcflags_);

  // Update R with obs errors that filter might have updated
  ObsVector_ obserr(Rmat_->obserrors());
  obserr = *obserrfilter_;
  Rmat_->update(obserr);

  // Save current obs, obs error estimates and QC flags (for diagnostics use only)
  std::string siter = "";
  if (iterconf_->has("iteration")) siter = iterconf_->getString("iteration");

  if (iterconf_->getBool("save qc", true)) {
    const std::string qcname = "EffectiveQC" + siter;
    qcflags_->save(qcname);
  }
  if (iterconf_->getBool("save hofx", true)) {
    const std::string obsname = "hofx" + siter;
    yobsim.save(obsname);
  }
  if (iterconf_->getBool("save obs errors", true)) {
    const std::string errname = "EffectiveError" + siter;
    obserrfilter_->save(errname);
  }
  if (iterconf_->getBool("save obs bias", true)) {
    const std::string biasname  = "ObsBias" + siter;
    ybias.save(biasname);
  }

  Log::info() << "Observer::finalize QC = " << *qcflags_ << std::endl;

  // Copy qc flags to pass out
  qcflags = *qcflags_;

  initialized_ = false;
  Log::trace() << "Observer<MODEL, OBS>::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void Observer<MODEL, OBS>::resetObsPert(const Geometry_ & geom,
                                        std::unique_ptr<ObsOperatorBase_> obsOpBase,
                                        const std::shared_ptr<GetValues_> & getValTL) {
  obsop_ = std::move(obsOpBase);
  getvals_ = getValTL;

  geovars_ = obsop_->requiredVars();
  varsizes_ = geom.variableSizes(geovars_);
  initialized_ = true;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
void Observer<MODEL, OBS>::updateObserver(const eckit::Configuration & cdaConfig) {
  getvals_->updateGetVals(cdaConfig);
  qcflags_->zeroAppended();
  obserrfilter_->readAppended("ObsError");
  Log::trace() << "Observer obs error appended" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVER_H_
