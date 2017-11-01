/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_POSTBASEAD_H_
#define OOPS_BASE_POSTBASEAD_H_

#include <boost/noncopyable.hpp>

#include "oops/base/PostTimer.h"
#include "eckit/config/Configuration.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Handles post-processing of model fields related to cost function.
/*!
 *  PostBaseAD is the base class for all cost function processors, it
 *  is mostly used so that PostProcessorAD can hold a vector of such
 *  processors.
 *  The difference with PostBase is that the adjoint can modify the
 *  increment.
 */

template <typename INCR> class PostBaseAD : private boost::noncopyable {
 public:
  PostBaseAD() : timer_() {}
  PostBaseAD(const util::DateTime & start, const util::DateTime & finish,
           const util::Duration & freq = util::Duration(0))
    : timer_(start, finish, freq) {}
  virtual ~PostBaseAD() {}

  void initializeAD(INCR & dx, const util::DateTime & bgn,
                    const util::Duration & step) {
    timer_.initialize(bgn, dx.validTime(), step);
    this->doFirstAD(dx, bgn, step);
  }

  void processAD(INCR & dx) {
    if (timer_.itIsTime(dx.validTime())) this->doProcessingAD(dx);
  }

  void finalizeAD(INCR & dx) {
    if (timer_.itIsTime(dx.validTime())) this->doProcessingAD(dx);
    this->doLastAD(dx);
  }

 private:
  PostTimer timer_;
  virtual void doFirstAD(INCR &, const util::DateTime &, const util::Duration &) =0;
  virtual void doProcessingAD(INCR &) =0;
  virtual void doLastAD(INCR &) =0;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_POSTBASEAD_H_
