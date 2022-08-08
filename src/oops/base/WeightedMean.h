/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_WEIGHTEDMEAN_H_
#define OOPS_BASE_WEIGHTEDMEAN_H_

#include <cmath>
#include <map>
#include <memory>

#include "eckit/exception/Exceptions.h"

#include "oops/base/Accumulator.h"
#include "oops/base/DolphChebyshev.h"
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

template <typename MODEL, typename FLDS>
class WeightedMean : public PostBase<FLDS> {
  typedef Geometry<MODEL>            Geometry_;

 public:
  WeightedMean(const Variables &, const util::DateTime &, const util::Duration &,
               const Geometry_ &, const eckit::Configuration &);
  virtual ~WeightedMean() {}

  FLDS * releaseMean();

 private:
  void doInitialize(const FLDS &, const util::DateTime &, const util::Duration &) override;

  void doProcessing(const FLDS &) override;

  std::unique_ptr<WeightingFct> wfct_;
  std::map< util::DateTime, double > weights_;
//  std::unique_ptr< Accumulator<MODEL, FLDS, FLDS> > avg_;
  Accumulator<MODEL, FLDS, FLDS> * avg_;
  double sum_;
  bool linit_;
  const util::DateTime bgn_;
  const util::DateTime end_;
  util::DateTime endleg_;
};

// =============================================================================

template <typename MODEL, typename FLDS>
WeightedMean<MODEL, FLDS>::WeightedMean(const Variables & vars,
                                        const util::DateTime & vt,
                                        const util::Duration & span,
                                        const Geometry_ & resol,
                                        const eckit::Configuration & config)
    : PostBase<FLDS>(vt-span/2, vt+span/2),
      wfct_(), weights_(), avg_(0), sum_(0.0), linit_(false),
      bgn_(vt-span/2), end_(vt+span/2), endleg_()
{
  avg_ = new Accumulator<MODEL, FLDS, FLDS>(resol, vars, vt);
  wfct_.reset(new DolphChebyshev(config));
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename FLDS>
FLDS * WeightedMean<MODEL, FLDS>::releaseMean() {
  ASSERT(linit_);
  ASSERT(std::abs(sum_ - 1.0) < 1.0e-8);
  return avg_;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename FLDS>
void WeightedMean<MODEL, FLDS>::doInitialize(const FLDS & xx,
                                             const util::DateTime & end,
                                             const util::Duration & tstep) {
  const util::DateTime bgn(xx.validTime());
  if (!linit_ && bgn <= end_ && end >= bgn_) {
    weights_ = wfct_->setWeights(bgn_, end_, tstep);
    linit_ = true;
  }
  endleg_ = end;
}

// -----------------------------------------------------------------------------

template <typename MODEL, typename FLDS>
void WeightedMean<MODEL, FLDS>::doProcessing(const FLDS & xx) {
  const util::DateTime now(xx.validTime());
  if (now != endleg_ || now == end_) {
    ASSERT(weights_.find(now) != weights_.end());
    const double zz = weights_[now];
    avg_->accumul(zz, xx);
    sum_ += zz;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_WEIGHTEDMEAN_H_
