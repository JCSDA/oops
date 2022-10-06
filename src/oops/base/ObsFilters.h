/*
 * (C) Copyright 2017-2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSFILTERS_H_
#define OOPS_BASE_OBSFILTERS_H_

#include <memory>
#include <set>
#include <vector>

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/ObsFilter.h"
#include "oops/base/ObsVector.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsSpace.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/parameters/Parameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/Printable.h"

namespace oops {

/// Configuration options for the ObsFilters class.
template <typename OBS>
class ObsFiltersParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsFiltersParameters, Parameters)

  typedef typename std::vector<oops::ObsFilterParametersWrapper<OBS>> FilterParams_;

 public:
  /// Options used to configure observation filters whose stage of operation
  /// (pre, prior or post) is determined automatically.
  oops::Parameter<FilterParams_> obsFilters{"obs filters", {}, this};

  /// Options used to configure observation pre filters.
  /// These filters are called before GetValues.
  /// Both GeoVaLs and H(x) are unavailable.
  oops::Parameter<FilterParams_> obsPreFilters{"obs pre filters", {}, this};

  /// Options used to configure observation filters.
  /// These filters are called after GetValues and before the observation operator.
  /// GeoVaLs are available and H(x) is unavailable.
  oops::Parameter<FilterParams_> obsPriorFilters{"obs prior filters", {}, this};

  /// Options used to configure observation filters.
  /// These filters are called after the observation operator.
  ///  Both GeoVaLs and H(x) are available.
  oops::Parameter<FilterParams_> obsPostFilters{"obs post filters", {}, this};
};

/// Holds observation filters (usually QC) for one observation type

// -----------------------------------------------------------------------------

template <typename OBS>
class ObsFilters : public util::Printable,
                   private boost::noncopyable {
  typedef GeoVaLs<OBS>            GeoVaLs_;
  typedef ObsDiagnostics<OBS>     ObsDiags_;
  typedef ObsSpace<OBS>           ObsSpace_;
  typedef ObsVector<OBS>          ObsVector_;
  typedef ObsFilter<OBS>          ObsFilter_;
  typedef ObsFiltersParameters<OBS> Parameters_;
  typedef std::shared_ptr<ObsDataVector<OBS, int> >  ObsDataPtr_;
  typedef std::vector<ObsFilterParametersWrapper<OBS>> FilterParams_;
  typedef ObsDataVector<OBS, float> ObsDataVector_;

 public:
  /// Initialize all filters for \p obspace, from parameters, using
  /// \p qcflags and \p obserr (observation error variances)
  /// \p iteration argument indicates outer loop iteration in the variational
  /// assimilation.
  ObsFilters(const ObsSpace_ &,
             const Parameters_ &,
             ObsDataPtr_ qcflags, ObsDataVector_ & obserr, const int iteration = 0);

  void preProcess();
  void priorFilter(const GeoVaLs_ &);
  void postFilter(const GeoVaLs_ &,
                  const ObsVector_ &,
                  const ObsVector_ &,
                  const ObsDiags_ &);

  Variables requiredVars() const {return geovars_;}
  Variables requiredHdiagnostics() const {return diagvars_;}

 private:
  void print(std::ostream &) const override;

  /// Configure a filter and append it to a vector of filters.
  void appendToFiltersList(const FilterParams_ & filtersParams,
                           std::vector<ObsFilter_> & filters);

  const ObsSpace_ & obsspace_;
  // List of filters for which the stage (pre/prior/post) will be determined automatically.
  std::vector<ObsFilter_> autoFilters_;
  // List of filters which have been designated to run at the pre stage.
  std::vector<ObsFilter_> preFilters_;
  // List of filters which have been designated to run at the prior stage.
  std::vector<ObsFilter_> priorFilters_;
  // List of filters which have been designated to run at the post stage .
  std::vector<ObsFilter_> postFilters_;
  Variables geovars_;
  Variables diagvars_;
  ObsDataPtr_ qcflags_;
  ObsDataVector_ & obserrfilter_;
  std::shared_ptr<ObsDataVector<OBS, float> > obserrtmp_;
  const int iteration_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsFilters<OBS>::ObsFilters(const ObsSpace_ & os,
                            const Parameters_ & params,
                            ObsDataPtr_ qcflags, ObsDataVector_ & obserr, const int iteration)
  : obsspace_(os), autoFilters_(), preFilters_(), priorFilters_(), postFilters_(),
    geovars_(), diagvars_(), qcflags_(qcflags), obserrfilter_(obserr),
    obserrtmp_(new ObsDataVector_(obserr)), iteration_(iteration) {
  Log::trace() << "ObsFilters::ObsFilters starting";

  const FilterParams_ & autoFiltersParams = params.obsFilters;
  const FilterParams_ & preFiltersParams = params.obsPreFilters;
  const FilterParams_ & priorFiltersParams = params.obsPriorFilters;
  const FilterParams_ & postFiltersParams = params.obsPostFilters;

  const bool atLeastOneAutoFilterConfigured = autoFiltersParams.size() > 0;
  const bool atLeastOnePrePriorPostFilterConfigured =
    preFiltersParams.size() +
    priorFiltersParams.size() +
    postFiltersParams.size() > 0;

  if (atLeastOneAutoFilterConfigured && atLeastOnePrePriorPostFilterConfigured) {
    throw eckit::UserError("It is not possible to use both an `obs filters` option "
                           "and one or more of the `obs pre filters`, "
                           "`obs prior filters` and `obs post filters` options", Here());
  }

  // If at least one filter has been configured, prepend the QC manager and
  // append the Final Check to the list of filters.
  const bool atLeastOneFilterConfigured =
    atLeastOneAutoFilterConfigured || atLeastOnePrePriorPostFilterConfigured;

// Prepare QC handling and statistics if any filters are present
  if (atLeastOneFilterConfigured) {
    eckit::LocalConfiguration conf;
    conf.set("filter", "QCmanager");
    ObsFilterParametersWrapper<OBS> filterParams;
    filterParams.deserialize(conf);
    if (atLeastOneAutoFilterConfigured)
      autoFilters_.emplace_back(os, filterParams.filterParameters, qcflags_, obserrtmp_);
    else
      preFilters_.emplace_back(os, filterParams.filterParameters, qcflags_, obserrtmp_);
  }

// Create the filters, only at 0-th iteration, or at iterations specified in "apply at iterations"
  if (atLeastOneAutoFilterConfigured)
    appendToFiltersList(autoFiltersParams, autoFilters_);
  if (atLeastOnePrePriorPostFilterConfigured) {
    appendToFiltersList(preFiltersParams, preFilters_);
    appendToFiltersList(priorFiltersParams, priorFilters_);
    appendToFiltersList(postFiltersParams, postFilters_);
  }

// Create the final filter run at the end of the pipeline
  if (atLeastOneFilterConfigured) {
    eckit::LocalConfiguration conf;
    conf.set("filter", "Final Check");
    ObsFilterParametersWrapper<OBS> filterParams;
    filterParams.deserialize(conf);
    if (atLeastOneAutoFilterConfigured)
      autoFilters_.emplace_back(os, filterParams.filterParameters, qcflags_, obserrtmp_);
    else
      postFilters_.emplace_back(os, filterParams.filterParameters, qcflags_, obserrtmp_);
  }

  Log::trace() << "ObsFilters::ObsFilters done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsFilters<OBS>::appendToFiltersList(const FilterParams_ & filtersParams,
                                          std::vector<ObsFilter_> & filters) {
  for (const ObsFilterParametersWrapper<OBS> & filterParams : filtersParams) {
    // Only create filters for the 0-th iteration by default
    bool apply = (iteration_ == 0);
    // If "apply at iterations" is set, check if this is the right iteration
    if (filterParams.applyAtIterations.value() != boost::none) {
      const std::set<int> iters = parseIntSet(*filterParams.applyAtIterations.value());
      apply = contains(iters, iteration_);
    }
    if (apply) {
      filters.emplace_back(obsspace_, filterParams.filterParameters, qcflags_, obserrtmp_);
      geovars_ += filters.back().requiredVars();
      diagvars_ += filters.back().requiredHdiagnostics();
    }
  }
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsFilters<OBS>::preProcess() {
  if (autoFilters_.size() > 0) {
    for (ObsFilter_ & filter : autoFilters_) {
      filter.checkFilterData(FilterStage::AUTO);
      filter.preProcess();
    }
  } else {
    for (ObsFilter_ & filter : preFilters_) {
      filter.checkFilterData(FilterStage::PRE);
      filter.preProcess();
    }
  }

  obserrfilter_ = *obserrtmp_;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsFilters<OBS>::priorFilter(const GeoVaLs_ & gv) {
  if (autoFilters_.size() > 0) {
    for (ObsFilter_ & filter : autoFilters_) {
      filter.checkFilterData(FilterStage::AUTO);
      filter.priorFilter(gv);
    }
  } else {
    for (ObsFilter_ & filter : priorFilters_) {
      filter.checkFilterData(FilterStage::PRIOR);
      filter.priorFilter(gv);
    }
  }

  obserrfilter_ = *obserrtmp_;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsFilters<OBS>::postFilter(const GeoVaLs_ & gv,
                                 const ObsVector_ & hofx,
                                 const ObsVector_ & bias,
                                 const ObsDiags_ & diags) {
  if (autoFilters_.size() > 0) {
    for (ObsFilter_ & filter : autoFilters_) {
      filter.checkFilterData(FilterStage::AUTO);
      filter.postFilter(gv, hofx, bias, diags);
    }
  } else {
    for (ObsFilter_ & filter : postFilters_) {
      filter.checkFilterData(FilterStage::POST);
      filter.postFilter(gv, hofx, bias, diags);
    }
  }

  obserrtmp_->mask(*qcflags_);
  obserrfilter_ = *obserrtmp_;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsFilters<OBS>::print(std::ostream & os) const {
  if (autoFilters_.size() > 0) {
    os << "ObsFilters: " << autoFilters_.size() << " elements:" << std::endl;
    for (const ObsFilter_ & filter : autoFilters_)
      os << filter << std::endl;
  } else {
    os << "Pre filters: " << preFilters_.size() << " elements:" << std::endl;
    for (const ObsFilter_ & preFilter : preFilters_)
      os << preFilter << std::endl;
    os << "Prior filters: " << priorFilters_.size() << " elements:" << std::endl;
    for (const ObsFilter_ & priorFilter : priorFilters_)
      os << priorFilter << std::endl;
    os << "Post filters: " << postFilters_.size() << " elements:" << std::endl;
    for (const ObsFilter_ & postFilter : postFilters_)
      os << postFilter << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSFILTERS_H_
