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
  : options_(), bgn_(), end_(), start_(), finish_() {}
// -----------------------------------------------------------------------------
PostTimer::PostTimer(const eckit::Configuration & conf)
  : options_(), bgn_(), end_(), start_(), finish_() {
  options_.deserialize(conf);
}
// -----------------------------------------------------------------------------
PostTimer::PostTimer(const util::DateTime & start, const util::DateTime & finish,
                     const util::Duration & freq)
  : options_(), bgn_(), end_(),
    start_(new util::DateTime(start)), finish_(new util::DateTime(finish)) {
  // setup config with passed frequency to init the options
  eckit::LocalConfiguration conf;
  conf.set("frequency", freq.toString());
  options_.deserialize(conf);
}
// -----------------------------------------------------------------------------
void PostTimer::initialize(const util::DateTime & bgn, const util::DateTime & end) {
  // TODO(someone): this can be simplified after weak-constraint subwindow
  // refactoring (start_ and finish_ can be removed; bgn_ and end_ changed to ptr)
  util::DateTime start(bgn);
  if (start_) {
    start = *start_;
  }
  bgn_ = start;
  util::DateTime finish(end);
  if (finish_) {
    finish = *finish_;
  }
  end_ = finish;
  // increase bgn_ value if needed
  bgn_ += options_.first;
}
// -----------------------------------------------------------------------------
bool PostTimer::itIsTime(const util::DateTime & now) {
  bool doit = false;

  const util::Duration & freq = options_.frequency;
  const std::vector<util::DateTime> & steps = options_.steps;

  if (now >= bgn_ && now <= end_) {
    // use at every step, and no prespecified steps?
    doit = (freq.toSeconds() == 0 && steps.empty());
    // frequency specified?
    if (!doit && freq.toSeconds() > 0) {
      const util::Duration dt = now - bgn_;
      doit = (dt >= util::Duration(0) && dt % freq == 0);
    }
    // steps are prespecified?
    if (!doit && !steps.empty()) {
      auto it = find(steps.begin(), steps.end(), now);
      doit = (it != steps.end());
    }
  }

  Log::trace() << "In PostTimer:itIsTime, time = " << now << ", doit = " << doit << std::endl;
  return doit;
}
// -----------------------------------------------------------------------------
}  // namespace oops
