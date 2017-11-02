/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#ifndef OOPS_BASE_FILTERBASE_H_
#define OOPS_BASE_FILTERBASE_H_

#include <boost/noncopyable.hpp>

#include "oops/interface/GeoVaLs.h"
#include "oops/interface/ObservationSpace.h"
#include "oops/interface/ObsVector.h"

namespace oops {

/// Base class for QC filters applied to observations

// -----------------------------------------------------------------------------

template <typename MODEL>
class FilterBase : private boost::noncopyable {
  typedef GeoVaLs<MODEL>             GeoVaLs_;
  typedef ObservationSpace<MODEL>    ObsSpace_;
  typedef ObsVector<MODEL>           ObsVector_;

 public:
  FilterBase() {}
  virtual ~FilterBase() {}

  virtual void preProcess(const GeoVaLs_ &, const ObsSpace_ &) const =0;
  virtual void postProcess(const GeoVaLs_ &, const ObsVector_ &, const ObsSpace_ &) const =0;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_FILTERBASE_H_
