/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_POSTBASETL_H_
#define OOPS_BASE_POSTBASETL_H_

#include <boost/noncopyable.hpp>

#include "oops/base/GeneralizedDepartures.h"
#include "oops/base/PostTimer.h"
#include "eckit/config/Configuration.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Handles post-processing of model fields related to cost function.
/*!
 *  PostBaseTL is the base class for all cost function processors, it
 *  is mostly used so that PostProcessorTL can hold a vector of such
 *  processors.
 *  By default processing is performed on every call.
 */

template <typename INCR> class PostBaseTL : private boost::noncopyable {
 public:
  PostBaseTL() : timer_() {}
  PostBaseTL(const util::DateTime & start, const util::DateTime & finish,
             const util::Duration & freq = util::Duration(0))
    : timer_(start, finish, freq) {}
  virtual ~PostBaseTL() {}

  void initializeTL(const INCR & dx, const util::DateTime & end,
                    const util::Duration & step) {
    timer_.initialize(dx.validTime(), end, step);
    this->doInitializeTL(dx, end, step);
    if (timer_.itIsTime(dx.validTime())) this->doProcessingTL(dx);
  }

  void processTL(const INCR & dx) {
    if (timer_.itIsTime(dx.validTime())) this->doProcessingTL(dx);
  }

  void finalizeTL(const INCR & dx) {
    this->doFinalizeTL(dx);
  }

/// Return TL dual space output
  virtual GeneralizedDepartures * releaseOutputFromTL() =0;

 private:
  PostTimer timer_;
  virtual void doInitializeTL(const INCR &, const util::DateTime &,
                              const util::Duration &) =0;
  virtual void doProcessingTL(const INCR &) =0;
  virtual void doFinalizeTL(const INCR &) =0;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_POSTBASETL_H_
