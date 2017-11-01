/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_WEIGHTEDDIFFTL_H_
#define OOPS_BASE_WEIGHTEDDIFFTL_H_

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <cmath>
#include <map>

#include "oops/base/Accumulator.h"
#include "oops/base/DolphChebyshev.h"
#include "oops/base/PostBaseTL.h"
#include "oops/base/WeightingFct.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Variables.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Compute time average of states or increments during linear model run.
/*!
 *  Derived classes will compute different types of averages (plain
 *  mean, various types of digital filters) by overwriting the weights
 *  computation method.
 *
 *  A lot of code here is common with WeightedDiff and even WeightedMean,
 *  the design could be improved to reduce code duplication. YT
 */

template <typename MODEL, typename INCR>
class WeightedDiffTL : public PostBaseTL<INCR> {
  typedef Geometry<MODEL>            Geometry_;
  typedef Variables<MODEL>           Variables_;

 public:
  WeightedDiffTL(const util::DateTime &, const util::Duration &,
                 const Geometry_ &, const eckit::Configuration &,
                 const util::Duration &, WeightingFct &);

  virtual ~WeightedDiffTL() {}

  INCR * releaseOutputFromTL() override;

 private:

  void doInitializeTL(const INCR &, const util::DateTime &, const util::Duration &) override;
  void doProcessingTL(const INCR &) override;
  void doFinalizeTL(const INCR &) override {}

  WeightingFct & wfct_;
  std::map< util::DateTime, double > weights_;
  Accumulator<MODEL, INCR, INCR> * avg_;  // Should be unique_ptr
  double sum_;
  bool linit_;
  const util::DateTime vtime_;
  const util::DateTime bgn_;
  const util::DateTime end_;
  util::Duration tstep_;
  util::DateTime bgnleg_;
  util::DateTime endleg_;
};

// =============================================================================

template <typename MODEL, typename INCR>
WeightedDiffTL<MODEL, INCR>::WeightedDiffTL(const util::DateTime & vt,
                                            const util::Duration & span,
                                            const Geometry_ & resol,
                                            const eckit::Configuration & config,
                                            const util::Duration & tstep,
                                            WeightingFct & wfct)
  : PostBaseTL<INCR>(vt-span/2, vt+span/2),
    wfct_(wfct), weights_(), avg_(0), sum_(0.0), linit_(false),
    vtime_(vt), bgn_(vt-span/2), end_(vt+span/2), tstep_(tstep),
    bgnleg_(), endleg_()
{
  Variables_ vars(config);
  avg_ = new Accumulator<MODEL, INCR, INCR>(resol, vars, vtime_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename INCR>
INCR * WeightedDiffTL<MODEL, INCR>::releaseOutputFromTL() {
  ASSERT(linit_);
  ASSERT(std::abs(sum_) < 1.0e-8);
  return avg_;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename INCR>
void WeightedDiffTL<MODEL, INCR>::doInitializeTL(const INCR & dx,
                                                 const util::DateTime & end,
                                                 const util::Duration & tstep) {
  const util::DateTime bgn(dx.validTime());
  ASSERT(bgn <= end);
  if (!linit_ && bgn <= end_ && end >= bgn_) {
    if (tstep_ == util::Duration(0)) tstep_ = tstep;
    ASSERT(tstep_ > util::Duration(0));
    weights_ = wfct_.setWeights(bgn_, end_, tstep_);
    linit_ = true;
    ASSERT(weights_.find(vtime_) != weights_.end());
    weights_[vtime_] -= 1.0;
  }
  bgnleg_ = bgn;
  endleg_ = end;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename INCR>
void WeightedDiffTL<MODEL, INCR>::doProcessingTL(const INCR & xx) {
  const util::DateTime now(xx.validTime());
  if (((bgnleg_ < end_ && endleg_ > bgn_) || bgnleg_ == endleg_) &&
      (now != endleg_ || now == end_ || now == bgnleg_)) {
    ASSERT(weights_.find(now) != weights_.end());
    const double zz = weights_[now];
    avg_->axpy(zz, xx, false);
    sum_ += zz;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_WEIGHTEDDIFFTL_H_
