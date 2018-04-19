/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_WEIGHTEDDIFFTLAD_H_
#define OOPS_BASE_WEIGHTEDDIFFTLAD_H_

#include <cmath>
#include <map>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "oops/base/Accumulator.h"
#include "oops/base/DolphChebyshev.h"
#include "oops/base/PostBaseTLAD.h"
#include "oops/base/Variables.h"
#include "oops/base/WeightingFct.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
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

template <typename MODEL>
class WeightedDiffTLAD : public PostBaseTLAD<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
  WeightedDiffTLAD(const util::DateTime &, const util::Duration &,
                 const Geometry_ &, const eckit::Configuration &,
                 const util::Duration &, WeightingFct &);
  WeightedDiffTLAD(const util::DateTime &, const util::Duration &,
                 const util::Duration &, WeightingFct &,
                 boost::shared_ptr<const Increment_>);

  virtual ~WeightedDiffTLAD() {}

  Increment_ * releaseOutputFromTL() override;

 private:
  void doInitializeTraj(const State_ &,
                        const util::DateTime &, const util::Duration &) override {}
  void doProcessingTraj(const State_ &) override {}
  void doFinalizeTraj(const State_ &) override {}

  void doInitializeTL(const Increment_ &,
                      const util::DateTime &, const util::Duration &) override;
  void doProcessingTL(const Increment_ &) override;
  void doFinalizeTL(const Increment_ &) override {}

  void doFirstAD(Increment_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessingAD(Increment_ &) override;
  void doLastAD(Increment_ &) override {}

  WeightingFct & wfct_;
  std::map< util::DateTime, double > weights_;
  boost::shared_ptr<const Increment_> forcing_;
  Accumulator<MODEL, Increment_, Increment_> * avg_;  // Should be unique_ptr
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

template <typename MODEL>
WeightedDiffTLAD<MODEL>::WeightedDiffTLAD(const util::DateTime & vt,
                                          const util::Duration & span,
                                          const Geometry_ & resol,
                                          const eckit::Configuration & config,
                                          const util::Duration & tstep,
                                          WeightingFct & wfct)
  : PostBaseTLAD<MODEL>(vt-span/2, vt+span/2),
    wfct_(wfct), weights_(), avg_(0), sum_(0.0), linit_(false),
    vtime_(vt), bgn_(vt-span/2), end_(vt+span/2), tstep_(tstep),
    bgnleg_(), endleg_()
{
  const Variables vars(config);
  avg_ = new Accumulator<MODEL, Increment_, Increment_>(resol, vars, vtime_);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
WeightedDiffTLAD<MODEL>::WeightedDiffTLAD(const util::DateTime & vt,
                                          const util::Duration & span,
                                          const util::Duration & tstep,
                                          WeightingFct & wfct,
                                          boost::shared_ptr<const Increment_> forcing)
  : PostBaseTLAD<MODEL>(vt-span/2, vt+span/2),
    wfct_(wfct), weights_(), forcing_(forcing), sum_(0.0), linit_(false),
    vtime_(vt), bgn_(vt-span/2), end_(vt+span/2), tstep_(tstep),
    bgnleg_(), endleg_()
{}

// -----------------------------------------------------------------------------

template <typename MODEL>
void WeightedDiffTLAD<MODEL>::doInitializeTL(const Increment_ & dx,
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

template <typename MODEL>
void WeightedDiffTLAD<MODEL>::doProcessingTL(const Increment_ & xx) {
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

template <typename MODEL>
Increment<MODEL> * WeightedDiffTLAD<MODEL>::releaseOutputFromTL() {
  ASSERT(linit_);
  ASSERT(std::abs(sum_) < 1.0e-8);
  return avg_;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void WeightedDiffTLAD<MODEL>::doFirstAD(Increment_ & dx,
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

template <typename MODEL>
void WeightedDiffTLAD<MODEL>::doProcessingAD(Increment_ & dx) {
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

#endif  // OOPS_BASE_WEIGHTEDDIFFTLAD_H_
