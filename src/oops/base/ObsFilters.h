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
  typedef ObsDataVector<MODEL, int> ObsVectorInt_;
  typedef ObservationSpace<MODEL>   ObsSpace_;
  typedef ObsVector<MODEL>          ObsVector_;

 public:
  ObsFilters(const ObsSpace_ &, const eckit::Configuration &, const Variables &);
  ObsFilters();
  ~ObsFilters();

  void priorFilter(const GeoVaLs_ &) const;
  void postFilter(const ObsVector_ &) const;

  const Variables & requiredGeoVaLs() const {return geovars_;}

 private:
  void print(std::ostream &) const;

  std::vector< boost::shared_ptr<ObsFilterBase_> > filters_;
  Variables geovars_;

  boost::shared_ptr<ObsVector_> obserr_;
  boost::shared_ptr<ObsVectorInt_> qcflags_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsFilters<MODEL>::ObsFilters(const ObsSpace_ & os, const eckit::Configuration & conf,
                              const Variables & observed)
  : filters_(), geovars_(), obserr_(), qcflags_() {
  Log::trace() << "ObsFilters::ObsFilters starting " << conf << std::endl;
  const int iter = conf.getInt("iteration");

// Initialize obs error values
  if (iter == 0) {
    ObsVector_ obserr(os, observed);
    obserr.read("ObsError");
    obserr.save("EffectiveError");
  }

// Get filters configuration
  std::vector<eckit::LocalConfiguration> confs;
  conf.get("ObsFilters", confs);

// Prepare storage for QC flags using PreQC filter (all work done in filter constructor)
  if (confs.size() > 0) {
    eckit::LocalConfiguration preconf;
    preconf.set("Filter", "PreQC");
    preconf.set("QCname", "EffectiveQC");
    preconf.set("observed", observed.variables());
    filters_.push_back(FilterFactory<MODEL>::create(os, preconf));
    obserr_.reset(new ObsVector_(os, observed));
    qcflags_.reset(new ObsVectorInt_(os, observed));
  }

// Create the filters
  for (std::size_t jj = 0; jj < confs.size(); ++jj) {
    bool apply = true;
    if (confs[jj].has("apply_at_iterations")) {
      std::set<int> iters = parseIntSet(confs[jj].getString("apply_at_iterations"));
      apply = contains(iters, iter);
    }
    if (apply) {
      confs[jj].set("QCname", "EffectiveQC");
      confs[jj].set("observed", observed.variables());
      boost::shared_ptr<ObsFilterBase_> tmp(FilterFactory<MODEL>::create(os, confs[jj]));
      geovars_ += tmp->requiredGeoVaLs();
      filters_.push_back(tmp);
    }
  }
  Log::trace() << "ObsFilters::ObsFilters done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsFilters<MODEL>::ObsFilters() : filters_(), geovars_(), obserr_(), qcflags_() {}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsFilters<MODEL>::~ObsFilters() {
  if (obserr_) {
    obserr_->read("EffectiveError");
    qcflags_->read("EffectiveQC");
    obserr_->mask(*qcflags_);
    obserr_->save("EffectiveError");
  }

  obserr_.reset();
  qcflags_.reset();
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
