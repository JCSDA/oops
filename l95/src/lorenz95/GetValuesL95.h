/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_GETVALUESL95_H_
#define LORENZ95_GETVALUESL95_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace lorenz95 {
  class GomL95;
  class LocsL95;
  class StateL95;
  class Resolution;

/// \brief used for getting state values at observation locations
///        (state at nearest gridpoint is used)
class GetValuesL95 : public util::Printable,
                     private util::ObjectCounter<GetValuesL95> {
 public:
  static const std::string classname() {return "lorenz95::GetValuesL95";}

  /// \brief computes indices resolidx_ of nearest gridpoints for all locations \p locs
  GetValuesL95(const Resolution &, const LocsL95 & locs);

  /// \brief fills in \p geovals for all observations in the timeframe (\p t1, \p t2],
  /// \p geovals are equal to the value of \p state at the nearest gridpoint
  void fillGeoVaLs(const StateL95 & state, const util::DateTime & t1,
                   const util::DateTime & t2, GomL95 & geovals) const;

 private:
  void print(std::ostream &) const;
  std::vector<int> resolidx_;          // indices in geometry for each of the locations
  std::vector<util::DateTime> times_;  // times of all the observations
};
// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_GETVALUESL95_H_
