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

#include <vector>

#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Handles timing of post-processing and similar actions
/*!
 *  By default processing is performed on every call.
 */

class PostTimer : private boost::noncopyable {
 public:
  PostTimer();
  explicit PostTimer(const util::Duration &);
  explicit PostTimer(const eckit::Configuration &);
  PostTimer(const util::DateTime &, const eckit::Configuration &);
  PostTimer(const util::DateTime &, const util::DateTime &, const util::Duration &);
  ~PostTimer() {}

  void initialize(const util::DateTime &, const util::DateTime &,
                  const util::Duration &);
  bool itIsTime(const util::DateTime &);

 private:
  const eckit::LocalConfiguration conf_;
  util::Duration freq_;
  util::DateTime bgn_;
  util::DateTime end_;
  boost::scoped_ptr<util::DateTime> start_;
  boost::scoped_ptr<util::DateTime> finish_;
  std::vector<util::DateTime> pptimes_;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_POSTTIMER_H_
