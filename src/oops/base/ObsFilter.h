/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_BASE_OBSFILTER_H_
#define OOPS_BASE_OBSFILTER_H_

#include <vector>

#include <boost/shared_ptr.hpp>

#include "oops/base/FilterBase.h"
#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsVector.h"
#include "util/Printable.h"

namespace oops {

/// Controls application of QC filters to observations

// -----------------------------------------------------------------------------

template<typename MODEL>
class ObsFilter : public util::Printable {
  typedef FilterBase<MODEL>          FilterBase_;
  typedef GeoVaLs<MODEL>             GeoVaLs_;
  typedef ObservationSpace<MODEL>    ObsSpace_;
  typedef ObsVector<MODEL>           ObsVector_;

 public:
  ObsFilter() : filters_(0) {}
  ObsFilter(const ObsFilter & pp): filters_(pp.filters_) {}
  ~ObsFilter() {}

  void enrollFilter(FilterBase_ *);
  void enrollFilter(boost::shared_ptr<FilterBase_>);

  void postFilter(const GeoVaLs_ &, const ObsVector_ &, const ObsSpace_ &) const;

 private:
  void print(std::ostream &) const;
  std::vector< boost::shared_ptr<FilterBase_> > filters_;
  ObsFilter operator= (const ObsFilter &);
};

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsFilter<MODEL>::enrollFilter(FilterBase_ * pp) {
  if (pp != 0) {
    boost::shared_ptr<FilterBase_> sp(pp);
    filters_.push_back(sp);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsFilter<MODEL>::enrollFilter(boost::shared_ptr<FilterBase_> pp) {
  if (pp != 0) filters_.push_back(pp);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void ObsFilter<MODEL>::postFilter(const GeoVaLs_ & gv, const ObsVector_ & ovec,
                                  const ObsSpace_ & obsdb) const {
  for (std::size_t jf = 0; jf < filters_.size(); ++jf) {
    filters_.at(jf)->postFilter(gv, ovec, obsdb);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void ObsFilter<MODEL>::print(std::ostream & os) const {
  os << "ObsFilter " << filters_.size() << " filters:" << std::endl;
  for (std::size_t jj = 0; jj < filters_.size(); ++jj) {
    os << *filters_[jj] << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_OBSFILTER_H_
