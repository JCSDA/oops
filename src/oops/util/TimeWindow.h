/*
 * (C) Crown copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_TIMEWINDOW_H_
#define OOPS_UTIL_TIMEWINDOW_H_

#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Printable.h"

namespace util {

/// \brief Indicates which of the time window bounds are treated as inclusive
/// when dealing with observation times.
enum class InclusiveWindowBound{LOWER, UPPER};

inline InclusiveWindowBound stringToWindowBound(const std::string & bound) {
  if (bound == "end") {
    return InclusiveWindowBound::UPPER;
  } else if (bound == "begin") {
    return InclusiveWindowBound::LOWER;
  } else {
    throw eckit::UserError("Invalid bound specified", Here());
  }
}

/// \brief This class stores the assimilation window upper and lower bounds and
/// provides several helper functions.
class TimeWindow : public util::Printable {
 public:
  explicit TimeWindow(const eckit::LocalConfiguration & conf);

  /// \brief Return a new TimeWindow object spanning `midpoint` +/- `halfwidth`.
  /// The code adjusts the bounds of the new time window in order to ensure they
  /// do not exceed those of the generating window.
  /// The new window has the same inclusive/exclusive bounds as the generating one.
  util::TimeWindow createSubWindow(const util::DateTime & midpoint,
                                   const util::Duration & halfwidth) const;

  /// \brief Return a new TimeWindow object spanning from `begin` to `end`.
  /// The code adjusts the bounds of the new time window in order to ensure they
  /// do not exceed those of the generating window.
  /// The new window has the same inclusive/exclusive bounds as the generating one.
  util::TimeWindow createSubWindow(const util::DateTime & begin,
                                   const util::DateTime & end) const;

  /// \brief Create a mask that indicates which entries in a vector of `util::DateTime` objects lie
  /// within the time window. A value of `true` indicates that an entry lies in the window.
  /// The choice of inclusive/exclusive window bounds are accounted for when producing this mask.
  std::vector<bool> createTimeMask(const std::vector<util::DateTime> &) const;

  /// \brief Create a mask that indicates which entries in a vector of `int64_t`,
  /// each of which represents a number of seconds relative to an epoch, lie within the time window.
  /// A value of `true` indicates that an entry lies in the window.
  /// The choice of inclusive/exclusive window bounds are accounted for when producing this mask.
  std::vector<bool> createTimeMask(const std::vector<int64_t> &) const;

  /// Return the start of the window, assuming an inclusive bound.
  /// The choice of inclusive/exclusive window bounds is disregarded when returning this value.
  util::DateTime start() const {return winbgn_;}

  /// Return the end of the window, assuming an inclusive bound.
  /// The choice of inclusive/exclusive window bounds is disregarded when returning this value.
  util::DateTime end() const {return winend_;}

  /// Return the length of the window, assuming inclusive lower and upper bounds.
  /// The choice of inclusive/exclusive window bounds is disregarded when returning this value.
  util::Duration length() const {return winend_ - winbgn_;}

  /// Return the midpoint of the window, assuming inclusive lower and upper bounds.
  /// The choice of inclusive/exclusive window bounds is disregarded when returning this value.
  util::DateTime midpoint() const {return winbgn_ + this->length() / 2;}

  /// \brief Set the epoch relative to which util::DateTime objects are converted to
  /// a number of seconds.
  void setEpoch(const util::DateTime &) const;

 private:
  /// Constructor for internal use only.
  TimeWindow(const util::DateTime & winbgn,
             const util::DateTime & winend,
             const InclusiveWindowBound inclusiveBound = InclusiveWindowBound::UPPER);

  /// \brief Provide information about the time window configuration.
  void print(std::ostream &) const;

  /// \brief Convert a vector of `util::DateTime` objects to a vector of `int64_t`
  /// objects, each of which represent the number of seconds relative to chosen epoch.
  std::vector<int64_t> convertDateTimesToEpochTimes
    (const std::vector<util::DateTime> &, const util::DateTime &) const;

  const util::DateTime winbgn_;
  const util::DateTime winend_;
  mutable util::DateTime epoch_;
  const InclusiveWindowBound incBound_;
};

}  // namespace util

#endif  // OOPS_UTIL_TIMEWINDOW_H_
