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
#include <memory>
#include <utility>

#include "oops/base/Accumulator.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/PostBaseTLAD.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/base/WeightedDiff.h"
#include "oops/base/WeightingFct.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

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
  WeightedDiffTLAD(const Variables &, const util::DateTime &, const util::Duration &,
                   const util::Duration &, const Geometry_ &, WeightingFct &);
  virtual ~WeightedDiffTLAD() {}

  Increment_ * releaseDiff() {return wdiff_.releaseDiff();}
  void setupTL(const Geometry_ &);
  void finalTL(Increment_ &);
  void setupAD(std::shared_ptr<const Increment_>);

 private:
  void doInitializeTraj(const State_ &,
                        const util::DateTime &, const util::Duration &) override;
  void doProcessingTraj(const State_ &) override;
  void doFinalizeTraj(const State_ &) override;

  void doInitializeTL(const Increment_ &,
                      const util::DateTime &, const util::Duration &) override;
  void doProcessingTL(const Increment_ &) override;
  void doFinalizeTL(const Increment_ &) override {}

  void doFirstAD(Increment_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessingAD(Increment_ &) override;
  void doLastAD(Increment_ &) override {}

  Variables vars_;
  WeightingFct & wfct_;
  WeightedDiff<MODEL, Increment_, State_> wdiff_;
  std::map< util::DateTime, double > weights_;
  std::shared_ptr<const Increment_> forcing_;
  std::unique_ptr<Accumulator<MODEL, Increment_, Increment_>> avg_;
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
WeightedDiffTLAD<MODEL>::WeightedDiffTLAD(const Variables & vars,
                                          const util::DateTime & vt,
                                          const util::Duration & span,
                                          const util::Duration & tstep,
                                          const Geometry_ & resol,
                                          WeightingFct & wfct)
  : PostBaseTLAD<MODEL>(vt-span/2, vt+span/2),
    vars_(vars), wfct_(wfct), wdiff_(vars, vt, span, tstep, resol, wfct_),
    weights_(), forcing_(), avg_(), sum_(0.0), linit_(false),
    vtime_(vt), bgn_(vt-span/2), end_(vt+span/2), tstep_(tstep),
    bgnleg_(), endleg_()
{
  Log::trace() << "WeightedDiffTLAD::WeightedDiffTLAD" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void WeightedDiffTLAD<MODEL>::doInitializeTraj(const State_ & xx,
                       const util::DateTime & end, const util::Duration & tstep) {
  Log::trace() << "WeightedDiffTLAD::doInitializeTraj start" << std::endl;
  wdiff_.initialize(xx, end, tstep);
  Log::trace() << "WeightedDiffTLAD::doInitializeTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void WeightedDiffTLAD<MODEL>::doProcessingTraj(const State_ & xx) {
  Log::trace() << "WeightedDiffTLAD::doProcessingTraj start" << std::endl;
  wdiff_.process(xx);
  Log::trace() << "WeightedDiffTLAD::doProcessingTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void WeightedDiffTLAD<MODEL>::doFinalizeTraj(const State_ & xx) {
  Log::trace() << "WeightedDiffTLAD::doFinalizeTraj start" << std::endl;
  wdiff_.finalize(xx);
  Log::trace() << "WeightedDiffTLAD::doFinalizeTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void WeightedDiffTLAD<MODEL>::setupTL(const Geometry_ & resol) {
  Log::trace() << "WeightedDiffTLAD::setupTL start" << std::endl;
  avg_.reset(new Accumulator<MODEL, Increment_, Increment_>(resol, vars_, vtime_));
  Log::trace() << "WeightedDiffTLAD::setupTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void WeightedDiffTLAD<MODEL>::doInitializeTL(const Increment_ & dx,
                                             const util::DateTime & end,
                                             const util::Duration & tstep) {
  Log::trace() << "WeightedDiffTLAD::doInitializeTL start" << std::endl;
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
  Log::trace() << "WeightedDiffTLAD::doInitializeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void WeightedDiffTLAD<MODEL>::doProcessingTL(const Increment_ & xx) {
  Log::trace() << "WeightedDiffTLAD::doProcessingTL start" << std::endl;
  const util::DateTime now(xx.validTime());
  if (((bgnleg_ < end_ && endleg_ > bgn_) || bgnleg_ == endleg_) &&
      (now != endleg_ || now == end_ || now == bgnleg_)) {
    ASSERT(weights_.find(now) != weights_.end());
    const double zz = weights_[now];
    avg_->axpy(zz, xx, false);
    sum_ += zz;
  }
  Log::trace() << "WeightedDiffTLAD::doProcessingTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void WeightedDiffTLAD<MODEL>::finalTL(Increment_ & out) {
  Log::trace() << "WeightedDiffTLAD::finalTL start" << std::endl;
  ASSERT(linit_);
  ASSERT(std::abs(sum_) < 1.0e-8);
  out = *avg_;
  Log::trace() << "WeightedDiffTLAD::finalTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void WeightedDiffTLAD<MODEL>::setupAD(std::shared_ptr<const Increment_> forcing) {
  Log::trace() << "WeightedDiffTLAD::setupAD start" << std::endl;
  forcing_ = forcing;
  Log::trace() << "WeightedDiffTLAD::setupAD done" << std::endl;
}

// -----------------------------------------------------------------------------


template <typename MODEL>
void WeightedDiffTLAD<MODEL>::doFirstAD(Increment_ & dx,
                                        const util::DateTime & bgn,
                                        const util::Duration & tstep) {
  Log::trace() << "WeightedDiffTLAD::doFirstAD start" << std::endl;
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
  Log::trace() << "WeightedDiffTLAD::doFirstAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void WeightedDiffTLAD<MODEL>::doProcessingAD(Increment_ & dx) {
  Log::trace() << "WeightedDiffTLAD::doProcessingAD start" << std::endl;
  ASSERT(forcing_);
  const util::DateTime now(dx.validTime());
  if (((bgnleg_ < end_ && endleg_ > bgn_) || bgnleg_ == endleg_) &&
      (now != endleg_ || now == end_ || now == bgnleg_)) {
    ASSERT(weights_.find(now) != weights_.end());
    const double zz = weights_[now];
    dx.axpy(zz, *forcing_, false);
    sum_ += zz;
  }
  Log::trace() << "WeightedDiffTLAD::doProcessingAD done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_WEIGHTEDDIFFTLAD_H_
