/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_POSTTIMER_H_
#define OOPS_BASE_POSTTIMER_H_

#include <memory>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/base/PostTimerParameters.h"
#include "oops/util/DateTime.h"

namespace eckit {
  class Configuration;
}

namespace util {
  class Duration;
}

namespace oops {

// -----------------------------------------------------------------------------

/// Handles timing of post-processing and similar actions
/*!
 *  By default processing is performed on every call.
 */

class PostTimer : private boost::noncopyable {
 public:
  PostTimer();
  explicit PostTimer(const eckit::Configuration &);
  explicit PostTimer(const PostTimerParameters &);
  PostTimer(const util::DateTime &, const util::DateTime &, const util::Duration &);

  void initialize(const util::DateTime &, const util::DateTime &);
  bool itIsTime(const util::DateTime &);

 private:
  util::DateTime bgn_;
  util::DateTime end_;
  std::unique_ptr<util::DateTime> start_;
  std::unique_ptr<util::DateTime> finish_;
  util::Duration frequency_;
  util::Duration first_;
  std::vector<util::Duration> steps_;
  std::vector<util::DateTime> times_;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_POSTTIMER_H_
