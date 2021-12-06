/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_WEIGHTEDDIFF_H_
#define OOPS_BASE_WEIGHTEDDIFF_H_

#include <cmath>
#include <map>

#include "oops/base/Accumulator.h"
#include "oops/base/Geometry.h"
#include "oops/base/PostBase.h"
#include "oops/base/Variables.h"
#include "oops/base/WeightingFct.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Compute time average of states or increments during model run.
/*!
 *  Derived classes will compute different types of averages (plain
 *  mean, various types of digital filters) by overwriting the weights
 *  computation method.
 */

template <typename MODEL, typename INCR, typename FLDS>
class WeightedDiff : public PostBase<FLDS> {
  typedef Geometry<MODEL>            Geometry_;

 public:
  WeightedDiff(const Variables &, const util::DateTime &, const util::Duration &,
               const util::Duration &, const Geometry_ &, WeightingFct &);
  virtual ~WeightedDiff() {}

  INCR * releaseDiff();

 private:
  void doInitialize(const FLDS &, const util::DateTime &, const util::Duration &) override;
  void doProcessing(const FLDS &) override;

  WeightingFct & wfct_;
  std::map< util::DateTime, double > weights_;
  Accumulator<MODEL, INCR, FLDS> * avg_;
  double sum_;
  bool linit_;
  const util::DateTime vtime_;
  const util::DateTime bgn_;
  const util::DateTime end_;
  util::Duration tstep_;
  util::DateTime bgnleg_;
  util::DateTime endleg_;
  util::DateTime current_;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename INCR, typename FLDS>
WeightedDiff<MODEL, INCR, FLDS>::WeightedDiff(const Variables & vars,
                                              const util::DateTime & vt,
                                              const util::Duration & span,
                                              const util::Duration & tstep,
                                              const Geometry_ & resol,
                                              WeightingFct & wfct)
  : PostBase<FLDS>(vt-span/2, vt+span/2),
    wfct_(wfct), weights_(), avg_(0), sum_(0.0), linit_(false),
    vtime_(vt), bgn_(vt-span/2), end_(vt+span/2), tstep_(tstep),
    bgnleg_(), endleg_(), current_()
{
  avg_ = new Accumulator<MODEL, INCR, FLDS>(resol, vars, vtime_);
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename INCR, typename FLDS>
INCR * WeightedDiff<MODEL, INCR, FLDS>::releaseDiff() {
  ASSERT(linit_);
  ASSERT(std::abs(sum_) < 1.0e-8);
  return avg_;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename INCR, typename FLDS>
void WeightedDiff<MODEL, INCR, FLDS>::doInitialize(const FLDS & xx,
                                                   const util::DateTime & end,
                                                   const util::Duration & tstep) {
  const util::DateTime bgn(xx.validTime());
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
  current_ = bgn-tstep;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename INCR, typename FLDS>
void WeightedDiff<MODEL, INCR, FLDS>::doProcessing(const FLDS & xx) {
  const util::DateTime now(xx.validTime());
  ASSERT(now > current_);
  if (((bgnleg_ < end_ && endleg_ > bgn_) || bgnleg_ == endleg_) &&
      (now != endleg_ || now == end_ || now == bgnleg_)) {
    ASSERT(weights_.find(now) != weights_.end());
    const double zz = weights_[now];
    avg_->accumul(zz, xx);
    sum_ += zz;
  }
  current_ = now;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_WEIGHTEDDIFF_H_
