// (C) Crown Copyright 2026 Met Office
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#include "oops/base/TimeInterpolationAccumulator.h"

#include "eckit/exception/Exceptions.h"

#include "oops/util/missingValues.h"

namespace oops {

namespace {

/// \brief Accumulate one state's time-weighted contribution into one observation column.
///
/// Both `rhs` and `lhs` address `nb_levs` contiguous values for a single observation.
///
/// \param isStateTime true if the observation time coincides with the state time, in
///                    which case the state's value is taken directly and no weighting
///                    is applied.
/// \param isFirst     true if this is the first of the (at most two) states contributing
///                    to this observation, i.e. the observation time falls after the
///                    state time. The first contribution assigns; the second adds.
/// \param timeWeight  fractional weight in [0, 1]; ignored when isStateTime.
/// \param nb_levs     number of vertical levels.
/// \param rhs         interpolated values for this observation [nb_levs].
/// \param lhs         accumulation target [nb_levs].
/// \param missing     the missingValue<double>() sentinel.
void timeAccumulateKernel(const bool isStateTime,
                          const bool isFirst,
                          const double timeWeight,
                          const size_t nb_levs,
                          const double * const rhs,
                          double * const lhs,
                          const double missing) {
  for (size_t jlev = 0; jlev < nb_levs; ++jlev) {
    if (rhs[jlev] == missing) {
      lhs[jlev] = missing;
    } else if (isStateTime) {
      lhs[jlev] = rhs[jlev];
    } else if (isFirst) {
      lhs[jlev] = rhs[jlev] * timeWeight;
    } else if (lhs[jlev] != missing) {
      // Don't accumulate onto a missing first contribution
      lhs[jlev] += rhs[jlev] * timeWeight;
    }
  }
}

}  // namespace

// -----------------------------------------------------------------------------

size_t TimeInterpolationAccumulator::accumulate(const std::vector<util::DateTime> & obsTimes,
                                                const std::vector<bool> & mask,
                                                const std::vector<size_t> & varLevels,
                                                const std::vector<double> & rhs,
                                                double * const lhs) const {
  ASSERT(mask.size() == obsTimes.size());

  const double missing = util::missingValue<double>();
  const size_t nb_obs = obsTimes.size();
  size_t offset = 0;

  for (const auto& nb_levs : varLevels) {
    for (size_t jloc = 0; jloc < nb_obs; ++jloc) {
      if (!mask[jloc]) continue;
      const util::DateTime & obsTime = obsTimes[jloc];
      timeAccumulateKernel(isStateTime(obsTime), isFirstContribution(obsTime),
                           weight(obsTime), nb_levs,
                           rhs.data() + offset + jloc * nb_levs,
                           lhs + offset + jloc * nb_levs,
                           missing);
    }
    offset += nb_obs * nb_levs;
  }

  ASSERT(offset == rhs.size());
  return offset;
}

}  // namespace oops
