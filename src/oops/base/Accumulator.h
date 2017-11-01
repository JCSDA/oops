/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_ACCUMULATOR_H_
#define OOPS_BASE_ACCUMULATOR_H_

#include "oops/interface/Geometry.h"
#include "oops/interface/Variables.h"
#include "util/DateTime.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL, typename ACC, typename FLDS> class Accumulator : public ACC {
  typedef Geometry<MODEL>            Geometry_;
  typedef Variables<MODEL>           Variables_;
 public:
  explicit Accumulator(const Geometry_ & resol, const Variables_ & vars,
                       const util::DateTime & vt)
    : ACC(resol, vars, vt) {ACC::zero();}
  ~Accumulator() {}

  void accumul(const double & zz, const FLDS & xx) {ACC::accumul(zz, xx);}
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_ACCUMULATOR_H_
