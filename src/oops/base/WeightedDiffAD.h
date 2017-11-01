/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_WEIGHTEDDIFFAD_H_
#define OOPS_BASE_WEIGHTEDDIFFAD_H_

#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <cmath>
#include <map>

#include "oops/base/DolphChebyshev.h"
#include "oops/base/PostBaseAD.h"
#include "oops/base/WeightingFct.h"
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

template <typename INCR>
class WeightedDiffAD : public PostBaseAD<INCR> {
 public:
  WeightedDiffAD(const util::DateTime &, const util::Duration &,
                 const util::Duration &, WeightingFct &,
                 boost::shared_ptr<const INCR>);

  virtual ~WeightedDiffAD() {}

 private:

  void doFirstAD(INCR &, const util::DateTime &, const util::Duration &) override;
  void doProcessingAD(INCR &) override;
  void doLastAD(INCR &) override {}

  WeightingFct & wfct_;
  std::map< util::DateTime, double > weights_;
  boost::shared_ptr<const INCR> forcing_;
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

template <typename INCR>
WeightedDiffAD<INCR>::WeightedDiffAD(const util::DateTime & vt,
                                     const util::Duration & span,
                                     const util::Duration & tstep,
                                     WeightingFct & wfct,
                                     boost::shared_ptr<const INCR> forcing)
  : PostBaseAD<INCR>(vt-span/2, vt+span/2),
    wfct_(wfct), weights_(), forcing_(forcing), sum_(0.0), linit_(false),
    vtime_(vt), bgn_(vt-span/2), end_(vt+span/2), tstep_(tstep),
    bgnleg_(), endleg_()
{}

// -----------------------------------------------------------------------------

template <typename INCR>
void WeightedDiffAD<INCR>::doFirstAD(INCR & dx,
                                     const util::DateTime & bgn,
                                     const util::Duration & tstep) {
  const util::DateTime end(dx.validTime());
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

template <typename INCR>
void WeightedDiffAD<INCR>::doProcessingAD(INCR & dx) {
  const util::DateTime now(dx.validTime());
  if (((bgnleg_ < end_ && endleg_ > bgn_) || bgnleg_ == endleg_) &&
      (now != endleg_ || now == end_ || now == bgnleg_)) {
    ASSERT(weights_.find(now) != weights_.end());
    const double zz = weights_[now];
    dx.axpy(zz, *forcing_, false);
    sum_ += zz;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_WEIGHTEDDIFFAD_H_
