// (C) Crown Copyright 2026 Met Office
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#pragma once

#include <cstddef>
#include <vector>

#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Accumulates the time-weighted contribution of one model state into an
///        observation buffer.
///
/// When observations are interpolated in time, each observation receives a contribution
/// from the two model states that bracket it, weighted by the observation's position
/// within the interval. This class computes the weight for a given observation
/// time and writes the weighted values into the destination buffer.
///
/// It deliberately knows nothing about MODEL, OBS or MPI, so the weighting and buffer
/// placement can be tested directly rather than only through a full model run.
class TimeInterpolationAccumulator {
 public:
  /// \param stateTime time of the model state being processed.
  /// \param halfWidth half-width of the time interval. With the state at the centre an
  ///                  observation exactly `halfWidth` away from `stateTime` receives
  ///                  zero weight.
  TimeInterpolationAccumulator(const util::DateTime & stateTime,
                               const util::Duration & halfWidth)
    : stateTime_(stateTime),
      tPrevious_(stateTime - halfWidth),
      tNext_(stateTime + halfWidth),
      dt_(static_cast<double>(halfWidth.toSeconds())) {}

  /// \brief True if the observation coincides with the state time, so the state's value
  ///        is used directly rather than blended.
  bool isStateTime(const util::DateTime & obsTime) const { return obsTime == stateTime_; }

  /// \brief True if this state is the first of the two to contribute to `obsTime`.
  ///
  /// States are processed in ascending time order, so the first state to contribute to
  /// an observation is the one that precedes it, i.e. the one for which the observation
  /// lies in the future.
  bool isFirstContribution(const util::DateTime & obsTime) const {
    return obsTime > stateTime_;
  }

  /// \brief Fractional weight in [0, 1] for `obsTime`; zero when `isStateTime`.
  double weight(const util::DateTime & obsTime) const {
    if (isStateTime(obsTime)) return 0.0;
    return isFirstContribution(obsTime)
      ? static_cast<double>((tNext_ - obsTime).toSeconds()) / dt_
      : static_cast<double>((obsTime - tPrevious_).toSeconds()) / dt_;
  }

  /// \brief Accumulate the interpolation result for one block of observations.
  ///
  /// `rhs` and the region of `lhs` written both use the block ordering
  ///   variable (slowest) -> observation (middle) -> level (fastest),
  /// covering `obsTimes.size()` observations for every variable in `varLevels`.
  /// Observations whose `mask` entry is false are left untouched in `lhs`.
  ///
  /// \param obsTimes  observation times for this block.
  /// \param mask      per-observation flag; false means the observation lies outside the
  ///                  time interval this state contributes to.
  /// \param varLevels varLevels[jf] = number of levels of variable `jf`.
  /// \param rhs       interpolation result for this block.
  /// \param lhs       start of the destination block; must hold at least `rhs.size()`
  ///                  values.
  /// \returns the number of values consumed from `rhs` and spanned in `lhs`.
  size_t accumulate(const std::vector<util::DateTime> & obsTimes,
                    const std::vector<bool> & mask,
                    const std::vector<size_t> & varLevels,
                    const std::vector<double> & rhs,
                    double * lhs) const;

 private:
  util::DateTime stateTime_;
  util::DateTime tPrevious_;
  util::DateTime tNext_;
  double dt_;
};

}  // namespace oops
