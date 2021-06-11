/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_GETVALUESTLAD_H_
#define LORENZ95_GETVALUESTLAD_H_

#include <ostream>
#include <string>
#include <vector>

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace lorenz95 {
  class GomL95;
  class IncrementL95;
  class LocsL95;
  class StateL95;
  class Resolution;

/// \brief used for getting state values at observation locations
///        and applying its TL & AD
class GetValuesTLAD : public util::Printable,
                     private util::ObjectCounter<GetValuesTLAD> {
 public:
  static const std::string classname() {return "lorenz95::GetValuesTLAD";}

  /// \brief computes indices resolidx_ of nearest gridpoints for all locations \p locs
  GetValuesTLAD(const Resolution &, const LocsL95 & locs, const eckit::Configuration &);

  /// \brief fills in \p geovals for all observations in the timeframe (\p t1, \p t2],
  /// \p geovals are equal to the value of \p state at the nearest gridpoint
  void setTrajectory(const StateL95 & state, const util::DateTime & t1,
                     const util::DateTime & t2, GomL95 & geovals);
  /// \brief fills in \p geovals for all observations in the timeframe (\p t1, \p t2],
  /// as a tangent-linear operator applied to \p inc
  void fillGeoVaLsTL(const IncrementL95 & inc, const util::DateTime & t1,
                     const util::DateTime & t2, GomL95 & geovals) const;
  /// \brief fills in \p inc as adjoint operator applied to \p geovals for all
  /// observations in the timeframe (\p t1, \p t2]
  void fillGeoVaLsAD(IncrementL95 & inc, const util::DateTime & t1,
                     const util::DateTime & t2, const GomL95 & geovals) const;

 private:
  void print(std::ostream &) const;
  std::vector<int> resolidx_;          // indices in geometry for each of the locations
  std::vector<util::DateTime> times_;  // times of all the observations
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_GETVALUESTLAD_H_
