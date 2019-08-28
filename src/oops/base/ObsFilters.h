/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSFILTERS_H_
#define OOPS_BASE_OBSFILTERS_H_

#include <set>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/ObsFilterBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/interface/ObsDiagnostics.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/IntSetParser.h"
#include "oops/util/Printable.h"

namespace oops {

/// Holds observation filters (usually QC) for one observation type

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsFilters : public util::Printable,
                   private boost::noncopyable {
  typedef GeoVaLs<MODEL>            GeoVaLs_;
  typedef ObsDiagnostics<MODEL>     ObsDiags_;
  typedef ObsFilterBase<MODEL>      ObsFilterBase_;
  typedef ObservationSpace<MODEL>   ObsSpace_;
  typedef ObsVector<MODEL>          ObsVector_;
  typedef boost::shared_ptr<ObsFilterBase<MODEL> >  ObsFilterPtr_;
  template <typename DATA> using ObsDataPtr_ = boost::shared_ptr<ObsDataVector<MODEL, DATA> >;

 public:
  ObsFilters(const ObsSpace_ &, const eckit::Configuration &,
             ObsDataPtr_<int> qcflags = ObsDataPtr_<int>(),
             ObsDataPtr_<float> obserr = ObsDataPtr_<float>());
  ObsFilters();
  ~ObsFilters();

  void preProcess() const;
  void priorFilter(const GeoVaLs_ &) const;
  void postFilter(const ObsVector_ &, const ObsDiags_ &) const;

  const Variables & requiredGeoVaLs() const {return geovars_;}
  const Variables & requiredHdiagnostics() const {return diagvars_;}

 private:
  void print(std::ostream &) const;

  std::vector<ObsFilterPtr_> filters_;
  Variables geovars_;
  Variables diagvars_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsFilters<MODEL>::ObsFilters(const ObsSpace_ & os, const eckit::Configuration & conf,
                              ObsDataPtr_<int> qcflags, ObsDataPtr_<float> obserr)
  : filters_(), geovars_(), diagvars_() {
  Log::trace() << "ObsFilters::ObsFilters starting " << conf << std::endl;

// Prepare QC handling and statistics if QC flags present
  if (qcflags) {
    eckit::LocalConfiguration preconf;
    preconf.set("Filter", "QCmanager");
    filters_.push_back(FilterFactory<MODEL>::create(os, preconf, qcflags, obserr));
  }

// Get filters configuration
  std::vector<eckit::LocalConfiguration> confs;
  conf.get("ObsFilters", confs);

// Create the filters
  for (std::size_t jj = 0; jj < confs.size(); ++jj) {
    bool apply = true;
    if (confs[jj].has("apply_at_iterations")) {
      std::set<int> iters = parseIntSet(confs[jj].getString("apply_at_iterations"));
      const int iter = conf.getInt("iteration");
      apply = contains(iters, iter);
    }
    if (apply) {
      ObsFilterPtr_ tmp(FilterFactory<MODEL>::create(os, confs[jj], qcflags, obserr));
      geovars_ += tmp->requiredGeoVaLs();
      diagvars_ += tmp->requiredHdiagnostics();
      filters_.push_back(tmp);
    }
  }

  Log::trace() << "ObsFilters::ObsFilters done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsFilters<MODEL>::ObsFilters() : filters_(), geovars_(), diagvars_() {}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsFilters<MODEL>::~ObsFilters() {
  Log::trace() << "ObsFilters::~ObsFilters destructed" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsFilters<MODEL>::preProcess() const {
  for (std::size_t jj = 0; jj < filters_.size(); ++jj) {
    filters_.at(jj)->preProcess();
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsFilters<MODEL>::priorFilter(const GeoVaLs_ & gv) const {
  for (std::size_t jj = 0; jj < filters_.size(); ++jj) {
    filters_.at(jj)->priorFilter(gv);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsFilters<MODEL>::postFilter(const ObsVector_ & hofx, const ObsDiags_ & diags) const {
  for (std::size_t jj = 0; jj < filters_.size(); ++jj) {
    filters_.at(jj)->postFilter(hofx, diags);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObsFilters<MODEL>::print(std::ostream & os) const {
  os << "ObsFilters: " << filters_.size() << " elements:" << std::endl;
  for (std::size_t jj = 0; jj < filters_.size(); ++jj) {
    os << *filters_.at(jj) << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSFILTERS_H_
