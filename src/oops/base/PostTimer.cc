/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "oops/base/PostTimer.h"

#include <algorithm>
#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {
// -----------------------------------------------------------------------------
PostTimer::PostTimer()
  : bgn_(), end_(), start_(), finish_(), frequency_(0), first_(0), steps_() {}
// -----------------------------------------------------------------------------
PostTimer::PostTimer(const PostTimerParameters & parameters)
  : bgn_(), end_(), start_(), finish_(),
    frequency_(parameters.frequency), first_(parameters.first), steps_(parameters.steps)
{}
// -----------------------------------------------------------------------------
PostTimer::PostTimer(const util::DateTime & start, const util::DateTime & finish,
                     const util::Duration & freq)
  : bgn_(), end_(),
    start_(new util::DateTime(start)), finish_(new util::DateTime(finish)),
    frequency_(freq), first_(0), steps_()
{}
// -----------------------------------------------------------------------------
void PostTimer::initialize(const util::DateTime & bgn, const util::DateTime & end) {
  bgn_ = bgn;
  end_ = end;
  if (start_) {
    bgn_ = *start_;
  }
  if (finish_) {
    end_ = *finish_;
  }
  // increase bgn_ value if needed
  bgn_ += first_;
}
// -----------------------------------------------------------------------------
bool PostTimer::itIsTime(const util::DateTime & now) {
  bool doit = false;

  if (now >= bgn_ && now <= end_) {
    // use at every step, and no prespecified steps?
    doit = (frequency_.toSeconds() == 0 && steps_.empty());
    // frequency specified?
    if (!doit && frequency_.toSeconds() > 0) {
      const util::Duration dt = now - bgn_;
      doit = (dt >= util::Duration(0) && dt % frequency_ == 0);
    }
    // steps are prespecified?
    if (!doit && !steps_.empty()) {
      auto it = find(steps_.begin(), steps_.end(), now);
      doit = (it != steps_.end());
    }
  }

  Log::trace() << "In PostTimer:itIsTime, time = " << now << ", doit = " << doit << std::endl;
  return doit;
}
// -----------------------------------------------------------------------------
}  // namespace oops
