/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "lorenz95/GetValuesTLAD.h"

#include <fstream>
#include <string>

#include "eckit/exception/Exceptions.h"
#include "oops/util/Logger.h"

#include "lorenz95/GomL95.h"
#include "lorenz95/IncrementL95.h"
#include "lorenz95/LocsL95.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/StateL95.h"

namespace lorenz95 {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
GetValuesTLAD::GetValuesTLAD(const Resolution & resol,
                             const LocsL95 & locs,
                             const eckit::Configuration & linearGetValuesConf)
  : resolidx_(locs.size()), times_(locs.times())
{
  // find indices of gridpoints nearest to all observations (resolidx_)
  const int npoints = resol.npoints();
  const double dres = static_cast<double>(npoints);
  for (size_t jobs = 0; jobs < locs.size(); ++jobs) {
    int ii = round(locs[jobs] * dres);
    ASSERT(ii >= 0 && ii <= npoints);
    if (ii == npoints) ii = 0;
    resolidx_[jobs] = ii;
  }
}
// -----------------------------------------------------------------------------
void GetValuesTLAD::setTrajectory(const StateL95 & state, const util::DateTime & t1,
                                  const util::DateTime & t2, GomL95 & vals) {
  const FieldL95 & field = state.getField();
  for (unsigned int jobs = 0; jobs < times_.size(); ++jobs) {
    // only fill in geovals for (t1, t2) timeslot
    if (times_[jobs] > t1 && times_[jobs] <= t2) {
      vals[jobs] = field[resolidx_[jobs]];
    }
  }
}
// -----------------------------------------------------------------------------
/// Get increment values at obs locations
// -----------------------------------------------------------------------------
void GetValuesTLAD::fillGeoVaLsTL(const IncrementL95 & inc, const util::DateTime & t1,
                                  const util::DateTime & t2, GomL95 & vals) const {
  const FieldL95 & field = inc.getField();
  for (unsigned int jobs = 0; jobs < times_.size(); ++jobs) {
    // only fill in geovals for (t1, t2) timeslot
    if (times_[jobs] > t1 && times_[jobs] <= t2) {
      vals[jobs] = field[resolidx_[jobs]];
    }
  }
}
// -----------------------------------------------------------------------------
void GetValuesTLAD::fillGeoVaLsAD(IncrementL95 & inc, const util::DateTime & t1,
                                  const util::DateTime & t2, const GomL95 & vals) const {
  FieldL95 & field = inc.getField();
  for (unsigned int jobs = 0; jobs < times_.size(); ++jobs) {
    // only fill in geovals for (t1, t2) timeslot
    if (times_[jobs] > t1 && times_[jobs] <= t2) {
      field[resolidx_[jobs]] += vals[jobs];
    }
  }
}
// -----------------------------------------------------------------------------
void GetValuesTLAD::print(std::ostream & os) const {
  os << "Nearest neighbor interpolation GetValues TL/AD ";
}
// -----------------------------------------------------------------------------

}  // namespace lorenz95
