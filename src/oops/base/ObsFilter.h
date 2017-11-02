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
#include "util/Logger.h"

namespace oops {

/// Controls application of QC filters to observations

template<typename MODEL>
class ObsFilter {
  typedef FilterBase<MODEL>          FilterBase_;
  typedef GeoVaLs<MODEL>             GeoVaLs_;
  typedef ObservationSpace<MODEL>    ObsSpace_;
  typedef ObsVector<MODEL>           ObsVector_;

 public:
  ObsFilter() : filters_(0) {}
  ObsFilter(const ObsFilter & pp): filters_(pp.filters_) {}
  ~ObsFilter() {}

  void enrollFilter(FilterBase_ * pp) {
    if (pp != 0) {
      boost::shared_ptr<FilterBase_> sp(pp);
      filters_.push_back(sp);
    }
  }

  void enrollFilter(boost::shared_ptr<FilterBase_> pp) {
    if (pp != 0) filters_.push_back(pp);
  }

  void preProcess(const GeoVaLs_ & gv, const ObsSpace_ & obsdb) const {
    for (std::size_t jf = 0; jf < filters_.size(); ++jf) {
      filters_.at(jf)->preProcess(gv, obsdb);
    }
  }

  void postProcess(const GeoVaLs_ & gv, const ObsVector_ & ovec,
                   const ObsSpace_ & obsdb) const {
    for (std::size_t jf = 0; jf < filters_.size(); ++jf) {
      filters_.at(jf)->postProcess(gv, ovec, obsdb);
    }
  }

 private:
  std::vector< boost::shared_ptr<FilterBase_> > filters_;
  ObsFilter operator= (const ObsFilter &);
};

}  // namespace oops

#endif  // OOPS_BASE_OBSFILTER_H_
