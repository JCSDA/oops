/*
 * (C) Crown copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>
#include <vector>

#include "oops/util/TimeWindow.h"

namespace util {

TimeWindow::TimeWindow(const eckit::LocalConfiguration & conf)
  : winbgn_(util::DateTime(conf.getString("begin"))),
    winend_(conf.has("length") ?
            winbgn_ + util::Duration(conf.getString("length")) :
            util::DateTime(conf.getString("end"))),
    epoch_(util::DateTime(1970, 1, 1, 0, 0, 0)),
    incBound_(stringToWindowBound(conf.getString("bound to include", "end")))
{
  if (conf.has("length") && conf.has("end")) {
    throw eckit::UserError("The time window configuration must define either "
                           "'length' or 'end', but not both.", Here());
  }
}

// Private constructor for internal use only.
TimeWindow::TimeWindow(const util::DateTime & winbgn,
                       const util::DateTime & winend,
                       const InclusiveWindowBound inclusiveBound)
  : winbgn_(winbgn),
    winend_(winend),
    epoch_(util::DateTime(1970, 1, 1, 0, 0, 0)),
    incBound_(inclusiveBound)
{}

util::TimeWindow TimeWindow::createSubWindow(const util::DateTime & begintime,
                                             const util::DateTime & endtime) const {
  const util::DateTime tLower = std::max(begintime, winbgn_);
  const util::DateTime tUpper = std::min(endtime, winend_);
  return util::TimeWindow(tLower, tUpper, incBound_);
}

util::TimeWindow TimeWindow::createSubWindow(const util::DateTime & midpoint,
                                             const util::Duration & halfwidth) const {
  return this->createSubWindow(midpoint - halfwidth, midpoint + halfwidth);
}

void TimeWindow::setEpoch(const util::DateTime & epoch) const {
  epoch_ = epoch;
}

std::vector<int64_t> TimeWindow::convertDateTimesToEpochTimes
(const std::vector<util::DateTime> & obsTimes, const util::DateTime & epoch) const {
  std::vector<int64_t> obsEpochTimes(obsTimes.size());
  for (size_t jobs = 0; jobs < obsTimes.size(); ++jobs) {
    obsEpochTimes[jobs] = (obsTimes[jobs] - epoch).toSeconds();
  }
  return obsEpochTimes;
}

std::vector<bool> TimeWindow::createTimeMask(const std::vector<util::DateTime> & obsTimes) const {
  const std::vector<int64_t> obsEpochTimes = this->convertDateTimesToEpochTimes(obsTimes, epoch_);
  return this->createTimeMask(obsEpochTimes);
}

std::vector<bool> TimeWindow::createTimeMask(const std::vector<int64_t> & obsTimes) const {
  const int64_t tBegin = (winbgn_ - epoch_).toSeconds();
  const int64_t tEnd = (winend_ - epoch_).toSeconds();

  std::vector<bool> timeMask(obsTimes.size());
  if (incBound_ == InclusiveWindowBound::LOWER) {
    for (size_t jobs = 0; jobs < obsTimes.size(); ++jobs) {
      timeMask[jobs] = obsTimes[jobs] >= tBegin && obsTimes[jobs] < tEnd;
    }
  } else {
    for (size_t jobs = 0; jobs < obsTimes.size(); ++jobs) {
      timeMask[jobs] = obsTimes[jobs] > tBegin && obsTimes[jobs] <= tEnd;
    }
  }
  return timeMask;
}

void TimeWindow::print(std::ostream & os) const {
  os << "TimeWindow: "
     << "start = " << winbgn_ << ", "
     << "end = " << winend_ << ", ";
  if (incBound_ == InclusiveWindowBound::LOWER) {
    os << "inclusive lower bound.";
  } else {
    os << "inclusive upper bound.";
  }
  os << std::endl;
}

}  // namespace util
