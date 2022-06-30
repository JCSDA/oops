/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_POSTBASETLAD_H_
#define OOPS_BASE_POSTBASETLAD_H_

#include <memory>

#include <boost/noncopyable.hpp>

#include "oops/base/Increment.h"
#include "oops/base/PostTimer.h"
#include "oops/base/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Handles post-processing of model fields related to cost function.
/*!
 *  PostBaseTLAD is the base class for all cost function linearized processors,
 *  it is mostly used so that PostProcessorTL can hold a vector of such
 *  processors.
 *  By default processing is performed on every call.
 *  An important difference with PostBase is that the adjoint can modify the
 *  increment.
 */

template <typename MODEL>
class PostBaseTLAD : private boost::noncopyable {
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
  PostBaseTLAD() : timer_() {}
  PostBaseTLAD(const util::DateTime & start, const util::DateTime & finish,
             const util::Duration & freq = util::Duration(0))
    : timer_(start, finish, freq) {}
  virtual ~PostBaseTLAD() {}

/// Set linearization state
  void initializeTraj(const State_ & xx, const util::DateTime & end,
                      const util::Duration & step) {
    timer_.initialize(xx.validTime(), end);
    this->doInitializeTraj(xx, end, step);
  }

  void processTraj(const State_ & xx) {
    if (timer_.itIsTime(xx.validTime())) this->doProcessingTraj(xx);
  }

  void finalizeTraj(const State_ & xx) {
    this->doFinalizeTraj(xx);
  }

/// Tangent linear methods
  void initializeTL(const Increment_ & dx, const util::DateTime & end,
                    const util::Duration & step) {
    timer_.initialize(dx.validTime(), end);
    this->doInitializeTL(dx, end, step);
  }

  void processTL(const Increment_ & dx) {
    if (timer_.itIsTime(dx.validTime())) this->doProcessingTL(dx);
  }

  void finalizeTL(const Increment_ & dx) {
    this->doFinalizeTL(dx);
  }

/// Adjoint methods
  void initializeAD(Increment_ & dx, const util::DateTime & bgn,
                    const util::Duration & step) {
    timer_.initialize(bgn, dx.validTime());
    this->doFirstAD(dx, bgn, step);
  }

  void processAD(Increment_ & dx) {
    if (timer_.itIsTime(dx.validTime())) this->doProcessingAD(dx);
  }

  void finalizeAD(Increment_ & dx) {
    this->doLastAD(dx);
  }

 private:
  PostTimer timer_;
  virtual void doInitializeTraj(const State_ &,
                                const util::DateTime &, const util::Duration &) = 0;
  virtual void doProcessingTraj(const State_ &) = 0;
  virtual void doFinalizeTraj(const State_ &) = 0;

  virtual void doInitializeTL(const Increment_ &,
                              const util::DateTime &, const util::Duration &) = 0;
  virtual void doProcessingTL(const Increment_ &) = 0;
  virtual void doFinalizeTL(const Increment_ &) = 0;

  virtual void doFirstAD(Increment_ &, const util::DateTime &, const util::Duration &) = 0;
  virtual void doProcessingAD(Increment_ &) = 0;
  virtual void doLastAD(Increment_ &) = 0;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_POSTBASETLAD_H_
