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
#include "oops/base/ObsError.h"
#include "oops/base/ObserverUtils.h"
#include "oops/base/ObsOperatorBase.h"
#include "oops/base/ObsVariables.h"
#include "oops/base/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsFilter.h"
#include "oops/interface/ObsOperator.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

/// Configuration options for the ObsFilters.
class ObsFiltersParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsFiltersParameters, Parameters)

 public:
  /// Options used to configure observation filters for toy models
  oops::OptionalParameter<eckit::LocalConfiguration> obsFilter{"obs filtering",
                  this};

  /// Options used to configure observation filters whose stage of operation
  /// (pre, prior or post) is determined automatically.
  oops::OptionalParameter<eckit::LocalConfiguration> obsFilters{"obs filters",
                  this};

  /// Options used to configure observation pre filters.
  /// These filters are called before GetValues.
  /// Both GeoVaLs and H(x) are unavailable.
  oops::OptionalParameter<eckit::LocalConfiguration> obsPreFilters{"obs pre filters",
                  this};

  /// Options used to configure observation filters.
  /// These filters are called after GetValues and before the observation operator.
  /// GeoVaLs are available and H(x) is unavailable.
  oops::OptionalParameter<eckit::LocalConfiguration> obsPriorFilters{"obs prior filters",
                  this};

  /// Options used to configure observation filters.
  /// These filters are called after the observation operator.
  ///  Both GeoVaLs and H(x) are available.
  oops::OptionalParameter<eckit::LocalConfiguration> obsPostFilters{"obs post filters",
                  this};
};

template <typename OBS>
class ObserverParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(ObserverParameters, Parameters)

 public:
  oops::RequiredParameter<eckit::LocalConfiguration> obsOperator{"obs operator", this};
  // Options used to configure filter.
  ObsFiltersParameters filtersParameters{this};

  oops::Parameter<eckit::LocalConfiguration> getValues{
    "get values", eckit::LocalConfiguration(), this};

  // Options used by ObserverTLAD. In the current design there is some overlap between the options
  // used by Observer and ObserverTLAD, so for now we include these options in ObserverParameters
  // to simplify the transition to Parameters. Ultimately, it will likely be a cleaner design to
  // separate out the options into ObserverParameters and ObserverTLADParameters.
  oops::Parameter<bool> monitoringOnly{"monitoring only", false, this};
  oops::OptionalParameter<eckit::LocalConfiguration> linearObsOperator{"linear obs operator", this};
};

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
  typedef ObserverParameters<OBS>      Parameters_;
  typedef ObsFilter<OBS>               ObsFilter_;
  typedef ObsOperator<OBS>             ObsOperator_;
  typedef ObsOperatorBase<OBS>         ObsOperatorBase_;
  typedef ObsSpace<OBS>                ObsSpace_;
  typedef ObsVector<OBS>               ObsVector_;
  typedef ObsDataVector<OBS, float>    ObsDataVector_;

 public:
/// \brief Initializes ObsOperators, Locations, and QC data
  Observer(const ObsSpace_ & obspace, const Parameters_ & params,
           std::unique_ptr<ObsOperatorBase_> obsOpBase = nullptr);

/// \brief Initializes variables, obs bias, obs filter (could be different for
/// different iterations
  std::vector<std::shared_ptr<GetValues_>> initialize(const Geometry_ &, const ObsAuxCtrl_ &,
                                                      ObsError_ &, const eckit::Configuration &);

/// \brief Computes H(x) from the filled in GeoVaLs
  void finalize(ObsVector_ &, ObsDataInt_ &);

  int resetObsPert(const Geometry_ &, std::unique_ptr<ObsOperatorBase_>,
                   const std::shared_ptr<GetValueTLADs_> &, int &);

 private:
  typedef std::vector<size_t> VariableSizes;

  Parameters_                       parameters_;
  const ObsSpace_ &                 obspace_;       // ObsSpace used in H(x)
  Variables                         allVars_;       // All required variables
  VariableSizes                     allVarSizes_;   // Sizes of these variables
  std::unique_ptr<ObsOperatorBase_> obsop_;         // Obs operator
  std::unique_ptr<Locations_>       locations_;     // Obs locations
  const ObsAuxCtrl_ *               biascoeff_;     // bias coefficients
  ObsError_ *                       Rmat_;          // Obs error covariance
  std::unique_ptr<ObsFilter_>       filter_;        // QC filter
  std::shared_ptr<ObsDataVector_>   obserrfilter_;  // Obs error std dev for processed variables
  // Instances of GetValues. Each receives a list of model variables and a set of paths along which
  // these variables should be interpolated. The interpolated values are stored in a single GeoVaLs
  // object (shared between all instances of GetValues).
  std::vector<std::shared_ptr<GetValues_>> getvals_;
  std::shared_ptr<ObsDataInt_>      qcflags_;       // QC flags (should not be a pointer)
  bool                              initialized_;
  std::unique_ptr<eckit::LocalConfiguration> iterconf_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename OBS>
Observer<MODEL, OBS>::Observer(const ObsSpace_ & obspace, const Parameters_ & params,
                               std::unique_ptr<ObsOperatorBase_> obsOpBase)
  : parameters_(params), obspace_(obspace), obsop_(),
    biascoeff_(nullptr), filter_(), qcflags_(), initialized_(false)
{
  Log::trace() << "Observer::Observer start" << std::endl;
  /// Set up observation operators
  if (obsOpBase == nullptr) {
    obsop_.reset(new ObsOperator_(obspace_, parameters_.obsOperator));
  } else {
    obsop_ = std::move(obsOpBase);
  }
  qcflags_.reset(new ObsDataInt_(obspace_, obspace_.obsvariables()));
  obserrfilter_.reset(new ObsDataVector_(obspace_, obspace_.obsvariables(), "ObsError"));
  Log::trace() << "Observer::Observer done" << std::endl;
}

// -----------------------------------------------------------------------------
template <typename MODEL, typename OBS>
std::vector<std::shared_ptr<GetValues<MODEL, OBS>>>
Observer<MODEL, OBS>::initialize(const Geometry_ & geom, const ObsAuxCtrl_ & biascoeff,
                                 ObsError_ & R, const eckit::Configuration & conf) {
  Log::trace() << "Observer<MODEL, OBS>::initialize start" << std::endl;
// Save information for finalize
  iterconf_.reset(new eckit::LocalConfiguration(conf));
  biascoeff_ = &biascoeff;
  Rmat_ = &R;

  // Set up QC filter and run preprocess
  const int iterfilt = iterconf_->getInt("iteration", 0);
  filter_.reset(new ObsFilter_(obspace_,
                               parameters_.filtersParameters.toConfiguration(),
                               qcflags_, obserrfilter_, iterfilt));
  filter_->preProcess();

  if (!initialized_) {
// Get the list of required variables
    oops::Variables geovars = obsop_->requiredVars();
    geovars += biascoeff_->requiredVars();
    geovars += filter_->requiredVars();

// Get the observation locations and their discretizations
    locations_ = std::make_unique<Locations_>(obsop_->locations());

  // Required variables grouped by the set of paths along which they'll be interpolated
    const std::vector<oops::Variables> groupedVars = groupVariablesByLocationSamplingMethod(
          geovars, *locations_);

// Get variable sizes (i.e. the numbers of levels in the corresponding GeoVaLs)
    const std::vector<VariableSizes> groupedVarSizes = variableSizes(groupedVars, geom);

    std::tie(allVars_, allVarSizes_) = mergeVariablesAndSizes(groupedVars, groupedVarSizes);

// Set up GetValues
    getvals_ = makeGetValuesVector(parameters_.getValues, geom,
                                   obspace_.timeWindow(),
                                   *locations_, groupedVars);
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
  GeoVaLs_ geovals = makeAndFillGeoVaLs(*locations_, allVars_, allVarSizes_, getvals_);

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
int Observer<MODEL, OBS>::resetObsPert(const Geometry_ & geom,
                                       std::unique_ptr<ObsOperatorBase_> obsOpBase,
                                       const std::shared_ptr<GetValueTLADs_> & getValTLs,
                                       int & index) {
  obsop_ = std::move(obsOpBase);
  int size = getvals_.size();
  for (int ii = 0; ii < size; ++ii) {
    getvals_[ii] = (*getValTLs)[ii + index];
  }

  allVars_ = obsop_->requiredVars();
  allVarSizes_ = geom.variableSizes(allVars_);
  initialized_ = true;

  return size + index;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSERVER_H_
