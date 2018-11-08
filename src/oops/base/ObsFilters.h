/*
 * (C) Copyright 2017-2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_OBSFILTERS_H_
#define OOPS_BASE_OBSFILTERS_H_

#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/ObsFilterBase.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/util/Printable.h"

namespace oops {

/// Holds observation filters (usually QC) for one observation type

// -----------------------------------------------------------------------------

template <typename MODEL>
class ObsFilters : public util::Printable {
  typedef GeoVaLs<MODEL>            GeoVaLs_;
  typedef ObsFilterBase<MODEL>      ObsFilterBase_;
  typedef ObservationSpace<MODEL>   ObsSpace_;
  typedef ObsVector<MODEL>          ObsVector_;

 public:
  ObsFilters(const ObsSpace_ &, const eckit::Configuration &);
  ObsFilters();
  ObsFilters(const ObsFilters &);
  ~ObsFilters();

  void priorFilter(const GeoVaLs_ &) const;
  void postFilter(const ObsVector_ &) const;

 private:
  void print(std::ostream &) const;
  std::vector< boost::shared_ptr<ObsFilterBase_> > filters_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsFilters<MODEL>::ObsFilters(const ObsSpace_ & os, const eckit::Configuration & conf) {
  std::vector<eckit::LocalConfiguration> confs;
  conf.get("ObsFilters", confs);
  filters_.resize(confs.size());
  for (std::size_t jj = 0; jj < confs.size(); ++jj) {
    filters_[jj].reset(FilterFactory<MODEL>::create(os, confs[jj]));
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsFilters<MODEL>::ObsFilters(): filters_() {}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsFilters<MODEL>::ObsFilters(const ObsFilters & other): filters_(other.filters_) {}

// -----------------------------------------------------------------------------

template <typename MODEL>
ObsFilters<MODEL>::~ObsFilters() {}

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
    os << *filters_[jj] << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSFILTERS_H_
