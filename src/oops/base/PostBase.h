/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_POSTBASE_H_
#define OOPS_BASE_POSTBASE_H_

#include <boost/noncopyable.hpp>

#include "oops/base/PostTimer.h"
#include "eckit/config/Configuration.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Handles post-processing of model fields.
/*!
 *  PostBase is the base class for all state post processors, it
 *  is mostly used so that PostProcessor can hold a vector of such
 *  processors.
 *  By default processing is performed on every call.
 */

template <typename FLDS> class PostBase : private boost::noncopyable {
 public:
/// Constructors and basic operators
  PostBase() : timer_() {}
  explicit PostBase(const util::Duration & freq) : timer_(freq) {}
  explicit PostBase(const eckit::Configuration & conf) : timer_(conf) {}
  PostBase(const util::DateTime & start, const eckit::Configuration & conf)
    : timer_(start, conf) {}
  PostBase(const util::DateTime & start, const util::DateTime & finish,
           const util::Duration & freq = util::Duration(0))
    : timer_(start, finish, freq) {}

  virtual ~PostBase() {}

/// Setup
  void initialize(const FLDS & xx, const util::DateTime & end,
                  const util::Duration & tstep) {
    timer_.initialize(xx.validTime(), end, tstep);
    this->doInitialize(xx, end, tstep);
    if (timer_.itIsTime(xx.validTime())) this->doProcessing(xx);
  }

/// Process state or increment
  void process(const FLDS & xx) {
    if (timer_.itIsTime(xx.validTime())) this->doProcessing(xx);
  }

/// Final
  void finalize(const FLDS & xx) {
    this->doFinalize(xx);
  }

 private:
  PostTimer timer_;

/// Actual processing
  virtual void doProcessing(const FLDS &) =0;
  virtual void doInitialize(const FLDS &, const util::DateTime &,
                            const util::Duration &) {}
  virtual void doFinalize(const FLDS &) {}
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_POSTBASE_H_
