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
#include <string>
#include <vector>

#include <boost/tokenizer.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "util/Logger.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/abor1_cpp.h"

namespace oops {
// -----------------------------------------------------------------------------
PostTimer::PostTimer()
  : conf_(), freq_(0), bgn_(), end_(), start_(0), finish_(0), pptimes_() {}
// -----------------------------------------------------------------------------
PostTimer::PostTimer(const util::Duration & freq)
  : conf_(), freq_(freq), bgn_(), end_(), start_(0), finish_(0), pptimes_() {}
// -----------------------------------------------------------------------------
PostTimer::PostTimer(const eckit::Configuration & conf)
  : conf_(conf), freq_(0), bgn_(), end_(), start_(0), finish_(0), pptimes_() {}
// -----------------------------------------------------------------------------
PostTimer::PostTimer(const util::DateTime & start, const eckit::Configuration & conf)
  : conf_(conf), freq_(0), bgn_(), end_(), start_(new util::DateTime(start)),
    finish_(0), pptimes_()
{}
// -----------------------------------------------------------------------------
PostTimer::PostTimer(const util::DateTime & start, const util::DateTime & finish,
                     const util::Duration & freq)
  : conf_(), freq_(freq), bgn_(), end_(), start_(new util::DateTime(start)),
    finish_(new util::DateTime(finish)), pptimes_() {}
// -----------------------------------------------------------------------------
void PostTimer::initialize(const util::DateTime & bgn, const util::DateTime & end,
                           const util::Duration &) {
  util::DateTime start(bgn);
  if (start_ != 0) {
    start = *start_;
  }
  bgn_ = start;
  util::DateTime finish(end);
  if (finish_ != 0) {
    finish = *finish_;
  }
  end_ = finish;

// User specified configuration
  Log::debug() << "PProc configuration is:"  << conf_ << std::endl;
// User specified frequency and start
  if (conf_.has("frequency"))
    freq_ = util::Duration(conf_.getString("frequency"));
  if (conf_.has("first")) {
    const util::Duration first(conf_.getString("first"));
    bgn_ += first;
  }
// User specified steps
  if (conf_.has("steps")) {
    const std::string steps = conf_.getString("steps");
    boost::tokenizer<> tok(steps);
    for (boost::tokenizer<>::iterator it = tok.begin(); it != tok.end(); ++it) {
      const util::Duration step(*it);
      const util::DateTime tt = start+step;
      if (bgn <= tt && tt <= end) pptimes_.push_back(tt);
    }
  }
}
// -----------------------------------------------------------------------------
bool PostTimer::itIsTime(const util::DateTime & now) {
  bool doit = false;
  if (now >= bgn_ && now <= end_) {
    doit = (freq_.toSeconds() == 0 && pptimes_.empty());
    if (!doit && freq_.toSeconds()>0) {
      const util::Duration dt = now - bgn_;
      doit = (dt >= util::Duration(0) && dt % freq_ == 0);
    }
    if (!doit && !pptimes_.empty()) {
      std::vector<util::DateTime>::iterator it;
      it = find(pptimes_.begin(), pptimes_.end(), now);
      doit = (it != pptimes_.end());
    }
  }
  Log::trace() << "In PostTimer:itIsTime, time = " << now << ", doit = " << doit << std::endl;
  return doit;
}
// -----------------------------------------------------------------------------
}  // namespace oops
