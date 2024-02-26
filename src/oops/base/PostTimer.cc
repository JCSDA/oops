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
  : bgn_(), end_(), start_(), finish_(), frequency_(0), first_(0), steps_(), times_() {}
// -----------------------------------------------------------------------------
PostTimer::PostTimer(const eckit::Configuration & config)
  : bgn_(), end_(), start_(), finish_(),
    frequency_(config.getString("frequency", "PT0S")),
    first_(config.getString("first", "PT0S")), steps_(), times_()
{
  if (config.has("times")) {
    std::vector<std::string> stimes = config.getStringVector("times");
    for (size_t js = 0; js < stimes.size(); ++js) {
      times_.push_back(util::DateTime(stimes[js]));
    }
  }
  if (config.has("steps")) {
    std::vector<std::string> ssteps = config.getStringVector("steps");
    for (size_t js = 0; js < ssteps.size(); ++js) {
      steps_.push_back(util::Duration(ssteps[js]));
    }
  }
}
// -----------------------------------------------------------------------------
PostTimer::PostTimer(const PostTimerParameters & parameters)
  : PostTimer(parameters.toConfiguration())
{}
// -----------------------------------------------------------------------------
PostTimer::PostTimer(const util::DateTime & start, const util::DateTime & finish,
                     const util::Duration & freq)
  : bgn_(), end_(),
    start_(new util::DateTime(start)), finish_(new util::DateTime(finish)),
    frequency_(freq), first_(0), steps_(), times_()
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
  for (const util::Duration & step : steps_) {
    times_.push_back(bgn + step);
  }
}
// -----------------------------------------------------------------------------
bool PostTimer::itIsTime(const util::DateTime & now) {
  bool doit = false;

  if (now >= bgn_ && now <= end_) {
    // use at every step, and no prespecified steps?
    doit = (frequency_.toSeconds() == 0 && times_.empty());
    // frequency specified?
    if (!doit && frequency_.toSeconds() > 0) {
      const util::Duration dt = now - bgn_;
      doit = (dt >= util::Duration(0) && dt % frequency_ == 0);
    }
    // steps are prespecified?
    if (!doit && !times_.empty()) {
      auto it = find(times_.begin(), times_.end(), now);
      doit = (it != times_.end());
    }
  }

  Log::trace() << "In PostTimer:itIsTime, time = " << now << ", doit = " << doit << std::endl;
  return doit;
}
// -----------------------------------------------------------------------------
}  // namespace oops
