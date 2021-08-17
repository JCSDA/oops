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
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObsSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Printable.h"

namespace oops {

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
  typedef std::shared_ptr<ObsDataVector<OBS, int> >  ObsDataPtr_;

 public:
  /// Initialize all filters for \p obspace, from parameters, using
  /// \p qcflags and \p obserr (observation error variances)
  /// \p iteration argument indicates outer loop iteration in the variational
  /// assimilation
  ObsFilters(const ObsSpace_ &, const std::vector<ObsFilterParametersWrapper<OBS>> &,
             ObsDataPtr_ qcflags, ObsVector_ & obserr, const int iteration = 0);

  void preProcess();
  void priorFilter(const GeoVaLs_ &);
  void postFilter(const ObsVector_ &, const ObsVector_ &, const ObsDiags_ &);

  Variables requiredVars() const {return geovars_;}
  Variables requiredHdiagnostics() const {return diagvars_;}

 private:
  void print(std::ostream &) const override;

  std::vector<ObsFilter_> filters_;
  Variables geovars_;
  Variables diagvars_;
  ObsDataPtr_ qcflags_;
  ObsVector_ & obserr_;
  std::shared_ptr<ObsDataVector<OBS, float> > obserrtmp_;
};

// -----------------------------------------------------------------------------

template <typename OBS>
ObsFilters<OBS>::ObsFilters(const ObsSpace_ & os,
                            const std::vector<ObsFilterParametersWrapper<OBS>> & filtersParams,
                            ObsDataPtr_ qcflags, ObsVector_ & obserr, const int iteration)
  : filters_(), geovars_(), diagvars_(), qcflags_(qcflags), obserr_(obserr),
    obserrtmp_(new ObsDataVector<OBS, float>(obserr)) {
  Log::trace() << "ObsFilters::ObsFilters starting:\n";
  for (const ObsFilterParametersWrapper<OBS> &filterParams : filtersParams)
    Log::trace() << "  " << filterParams << std::endl;

// Prepare QC handling and statistics if any filters are present
  if (filtersParams.size() > 0) {
    eckit::LocalConfiguration conf;
    conf.set("filter", "QCmanager");
    ObsFilterParametersWrapper<OBS> filterParams;
    filterParams.validateAndDeserialize(conf);
    filters_.emplace_back(os, filterParams.filterParameters, qcflags_, obserrtmp_);
  }

// Create the filters, only at 0-th iteration, or at iterations specified in "apply at iterations"
  for (const ObsFilterParametersWrapper<OBS> &filterParams : filtersParams) {
    // Only create filters for the 0-th iteration by default
    bool apply = (iteration == 0);
    // If "apply at iterations" is set, check if this is the right iteration
    if (filterParams.applyAtIterations.value() != boost::none) {
      std::set<int> iters = parseIntSet(*filterParams.applyAtIterations.value());
      apply = contains(iters, iteration);
    }
    if (apply) {
      filters_.emplace_back(os, filterParams.filterParameters, qcflags_, obserrtmp_);
      geovars_ += filters_.back().requiredVars();
      diagvars_ += filters_.back().requiredHdiagnostics();
    }
  }

// Create the final filter run at the end of the pipeline
  if (filtersParams.size() > 0) {
    eckit::LocalConfiguration conf;
    conf.set("filter", "Final Check");
    ObsFilterParametersWrapper<OBS> filterParams;
    filterParams.validateAndDeserialize(conf);
    filters_.emplace_back(os, filterParams.filterParameters, qcflags_, obserrtmp_);
  }

  Log::trace() << "ObsFilters::ObsFilters done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsFilters<OBS>::preProcess() {
  for (ObsFilter_ & filter : filters_) {
    filter.preProcess();
  }
  obserrtmp_->mask(*qcflags_);
  obserr_ = *obserrtmp_;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsFilters<OBS>::priorFilter(const GeoVaLs_ & gv) {
  for (ObsFilter_ & filter : filters_) {
    filter.priorFilter(gv);
  }
  obserrtmp_->mask(*qcflags_);
  obserr_ = *obserrtmp_;
}

// -----------------------------------------------------------------------------

template<typename OBS>
void ObsFilters<OBS>::postFilter(const ObsVector_ & hofx, const ObsVector_ & bias,
                                 const ObsDiags_ & diags) {
  for (ObsFilter_ & filter : filters_) {
    filter.postFilter(hofx, bias, diags);
  }
  obserrtmp_->mask(*qcflags_);
  obserr_ = *obserrtmp_;
}

// -----------------------------------------------------------------------------

template <typename OBS>
void ObsFilters<OBS>::print(std::ostream & os) const {
  os << "ObsFilters: " << filters_.size() << " elements:" << std::endl;
  for (const ObsFilter_ & filter : filters_) {
    os << filter << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSFILTERS_H_
