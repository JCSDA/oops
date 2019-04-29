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
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/ObsFilterBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObsDataVector.h"
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
  typedef ObsFilterBase<MODEL>      ObsFilterBase_;
  typedef ObservationSpace<MODEL>   ObsSpace_;
  typedef ObsVector<MODEL>          ObsVector_;
  typedef boost::shared_ptr<ObsFilterBase<MODEL> >  ObsFilterPtr_;
  template <typename DATA> using ObsDataPtr_ = boost::shared_ptr<ObsDataVector<MODEL, DATA> >;

 public:
  ObsFilters(const ObsSpace_ &, const eckit::Configuration &, const Variables &,
             ObsDataPtr_<int> qcflags = ObsDataPtr_<int>(),
             ObsDataPtr_<float> obserr = ObsDataPtr_<float>());
  ObsFilters();
  ~ObsFilters();

  void priorFilter(const GeoVaLs_ &) const;
  void postFilter(const ObsVector_ &) const;

  const Variables & requiredGeoVaLs() const {return geovars_;}

 private:
  void print(std::ostream &) const;

  std::vector<ObsFilterPtr_> filters_;
  Variables geovars_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsFilters<MODEL>::ObsFilters(const ObsSpace_ & os, const eckit::Configuration & conf,
                              const Variables & observed,
                              ObsDataPtr_<int> qcflags, ObsDataPtr_<float> obserr)
  : filters_(), geovars_() {
  Log::trace() << "ObsFilters::ObsFilters starting " << conf << std::endl;

// Prepapre QC
  if (conf.getString("PreQC", "off") == "on") {
    eckit::LocalConfiguration preconf;
    preconf.set("Filter", "PreQC");
    preconf.set("observed", observed.variables());
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
      confs[jj].set("observed", observed.variables());
      ObsFilterPtr_ tmp(FilterFactory<MODEL>::create(os, confs[jj], qcflags, obserr));
      geovars_ += tmp->requiredGeoVaLs();
      filters_.push_back(tmp);
    }
  }

  Log::trace() << "ObsFilters::ObsFilters done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsFilters<MODEL>::ObsFilters() : filters_(), geovars_() {}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsFilters<MODEL>::~ObsFilters() {
  Log::trace() << "ObsFilters::~ObsFilters destructed" << std::endl;
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
void ObsFilters<MODEL>::postFilter(const ObsVector_ & ovec) const {
  for (std::size_t jj = 0; jj < filters_.size(); ++jj) {
    filters_.at(jj)->postFilter(ovec);
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
