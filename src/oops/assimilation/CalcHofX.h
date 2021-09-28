/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#ifndef OOPS_ASSIMILATION_CALCHOFX_H_
#define OOPS_ASSIMILATION_CALCHOFX_H_

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/ObsAuxControls.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsFilters.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/Locations.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsOperator.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

template <typename OBS>
class CalcHofXParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(CalcHofXParameters, Parameters)

  typedef typename OBS::ObsOperator::Parameters_ ObsOperatorParameters_;

 public:
  oops::RequiredParameter<ObsOperatorParameters_> obsOperator{"obs operator", this};
  oops::Parameter<std::vector<ObsFilterParametersWrapper<OBS>>> obsFilters{"obs filters", {}, this};
};

// -----------------------------------------------------------------------------

/// \brief Computes observation operator (from GeoVaLs), applies bias correction
///        and runs QC filters
template <typename OBS>
class CalcHofX {
  typedef GeoVaLs<OBS>               GeoVaLs_;
  typedef Locations<OBS>             Locations_;
  typedef ObsAuxControls<OBS>        ObsAuxCtrls_;
  typedef ObsDiagnostics<OBS>        ObsDiags_;
  typedef Observations<OBS>          Observations_;
  typedef ObsFilters<OBS>            ObsFilters_;
  typedef ObsOperator<OBS>           ObsOperator_;
  typedef ObsSpaces<OBS>             ObsSpaces_;
  typedef ObsVector<OBS>             ObsVector_;
  template <typename DATA> using ObsData_ = ObsDataVector<OBS, DATA>;
  template <typename DATA> using ObsDataVec_ = std::vector<std::shared_ptr<ObsData_<DATA>>>;

  typedef std::vector<std::unique_ptr<GeoVaLs_>>       GeoVaLsVec_;
  typedef std::vector<std::unique_ptr<Locations_>>     LocationsVec_;
  typedef std::vector<std::unique_ptr<ObsFilters_>>    ObsFiltersVec_;
  typedef std::vector<std::unique_ptr<ObsOperator_>>   ObsOperatorVec_;
  typedef std::vector<std::shared_ptr<ObsVector_>>     ObsVectorVec_;
  typedef std::vector<Variables>                       VariablesVec_;

 public:
/// \brief Initializes ObsOperators, Locations, and QC data
  CalcHofX(const ObsSpaces_ &, const eckit::Configuration &);

/// \brief Initializes variables, obs bias, obs filters (could be different for
/// different iterations
  void initialize(const ObsAuxCtrls_ &, const int iteration = 0);

/// \brief Computes H(x) from the filled in GeoVaLs
  Observations_ compute(const GeoVaLsVec_ &);

/// \brief accessor to the locations
  const LocationsVec_ & locations() const { return locations_; }
/// \brief accessor to variables required from the model
  const VariablesVec_ & requiredVars() const { return geovars_; }
/// \brief accessor to QC flags
  const ObsDataVec_<int> & qcflags() const { return qcflags_; }

/// \brief read QC flags from \p qcname variable from file
  void readQcFlags(const std::string & qcname);
/// \brief reset QC flags and ObsErrors
  void resetQc();
/// \brief save QC flags to ObsSpaces
  void saveQcFlags(const std::string &) const;
/// \brief save obs error variances (modified in QC) to ObsSpaces
  void saveObsErrors(const std::string &) const;
/// \brief mask obs errors with QC flags
  void maskObsErrors();

 private:
  eckit::LocalConfiguration obsconfig_;
  const ObsSpaces_ &   obspaces_;   // ObsSpaces used in H(x)
  ObsOperatorVec_      obsops_;     // Obs operators
  LocationsVec_        locations_;  // locations
  const ObsAuxCtrls_ * biascoeff_;  // bias coefficients
  ObsDataVec_<int>     qcflags_;    // QC flags
  ObsVectorVec_        obserrs_;    // Obs error variances (used in QC filters)
  ObsFiltersVec_       filters_;    // QC filters
  VariablesVec_        geovars_;    // variables required from the model
};

// -----------------------------------------------------------------------------

template <typename OBS>
CalcHofX<OBS>::CalcHofX(const ObsSpaces_ & obspaces,
                        const eckit::Configuration & config) :
  obsconfig_(config), obspaces_(obspaces), obsops_(), locations_(),
  biascoeff_(nullptr), qcflags_(), obserrs_(), filters_(),
  geovars_(obspaces_.size())
{
  std::vector<eckit::LocalConfiguration> obsconfs = config.getSubConfigurations();
  for (size_t jj = 0; jj < obspaces_.size(); ++jj) {
    // obsconfs[jj] contains not only options controlling the obs operator and filters (known to
    // CalcHofX) but also those controlling the obs space (unknown to it). So we can't call
    // validateAndDeserialize() here, since "obs space" would be treated as an unrecognized
    // keyword. In the long term the code constructing the CalcHofX will probably need to split
    // the contents of the "observations" vector into two vectors, one containing the "obs space"
    // sections and the other the "obs operator" and "obs filters" sections, and pass the former to
    // the constructor of ObsSpaces and the latter to the constructor of CalcHofX.
    CalcHofXParameters<OBS> observerParams;
    observerParams.deserialize(obsconfs[jj]);
    /// Set up observation operators
    obsops_.emplace_back(new ObsOperator_(obspaces_[jj], observerParams.obsOperator));
    locations_.emplace_back(new Locations_(obsops_[jj]->locations()));

    /// Allocate QC flags
    qcflags_.emplace_back(std::make_shared<ObsData_<int>>(obspaces[jj],
                            obspaces[jj].obsvariables()));
    /// Allocate and read initial obs error
    obserrs_.emplace_back(std::make_shared<ObsVector_>(obspaces[jj], "ObsError"));
  }
  Log::trace() << "CalcHofX<OBS> constructed" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void CalcHofX<OBS>::initialize(const ObsAuxCtrls_ & obsaux, const int iteration) {
  std::vector<eckit::LocalConfiguration> obsconfs = obsconfig_.getSubConfigurations();
  biascoeff_ = &obsaux;
  filters_.clear();
  for (size_t jj = 0; jj < obspaces_.size(); ++jj) {
    CalcHofXParameters<OBS> observerParams;
    observerParams.deserialize(obsconfs[jj]);
    /// Set up QC filters and run preprocess
    filters_.emplace_back(new ObsFilters_(obspaces_[jj], observerParams.obsFilters,
                                          qcflags_[jj], *obserrs_[jj], iteration));
    filters_[jj]->preProcess();

    /// Set up variables requested from the model
    geovars_[jj] += obsops_[jj]->requiredVars();
    geovars_[jj] += (*biascoeff_)[jj].requiredVars();
    geovars_[jj] += filters_[jj]->requiredVars();
  }
  Log::trace() << "CalcHofX<OBS>::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename OBS>
Observations<OBS> CalcHofX<OBS>::compute(const GeoVaLsVec_ & geovals) {
  oops::Log::trace() << "CalcHofX<OBS>::compute start" << std::endl;

  Observations<OBS> yobs(obspaces_);
  Observations<OBS> ybias(obspaces_);
  ybias.zero();
  for (size_t jj = 0; jj < obspaces_.size(); ++jj) {
    /// call prior filters
    filters_[jj]->priorFilter(*geovals[jj]);
    /// compute H(x)
    oops::Variables vars;
    vars += filters_[jj]->requiredHdiagnostics();
    vars += (*biascoeff_)[jj].requiredHdiagnostics();
    ObsDiags_ ydiags(obspaces_[jj], *locations_[jj], vars);
    obsops_[jj]->simulateObs(*geovals[jj], yobs[jj], (*biascoeff_)[jj], ybias[jj], ydiags);
    /// save bc
    ybias[jj].save("ObsBias");
    /// call posterior filters
    filters_[jj]->postFilter(yobs[jj], ybias[jj], ydiags);
  }
  oops::Log::trace() << "CalcHofX<OBS>::compute done" << std::endl;
  return yobs;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void CalcHofX<OBS>::readQcFlags(const std::string & qcname) {
  for (const auto & qcflag : qcflags_) {
    qcflag->read(qcname);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
void CalcHofX<OBS>::resetQc() {
  for (const auto & qcflag : qcflags_) {
    qcflag->zero();
  }
  for (const auto & obserr : obserrs_) {
    obserr->read("ObsError");
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
void CalcHofX<OBS>::saveQcFlags(const std::string & name) const {
  for (const auto & qcflag : qcflags_) {
    qcflag->save(name);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
void CalcHofX<OBS>::saveObsErrors(const std::string & name) const {
  for (const auto & obserr : obserrs_) {
    obserr->save(name);
  }
}

// -----------------------------------------------------------------------------

template <typename OBS>
void CalcHofX<OBS>::maskObsErrors() {
  for (size_t jj = 0; jj < obserrs_.size(); ++jj) {
    obserrs_[jj]->mask(*qcflags_[jj]);
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_CALCHOFX_H_
