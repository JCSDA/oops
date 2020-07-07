/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QG_MODEL_GETVALUESTLAD_H_
#define QG_MODEL_GETVALUESTLAD_H_

#include <ostream>
#include <string>

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "oops/qg/QgFortran.h"

namespace qg {
  class GeometryQG;
  class GomQG;
  class IncrementQG;
  class LocationsQG;
  class StateQG;

/// \brief used for getting state values at observation locations
//  and applying its TL & AD
class GetValuesTLAD : public util::Printable,
                      private util::ObjectCounter<GetValuesTLAD> {
 public:
  static const std::string classname() {return "qg::GetValuesTLAD";}

/// \brief saves all locations \p locs to use during filling GeoVaLs
  GetValuesTLAD(const GeometryQG &, const LocationsQG & locs);
  virtual ~GetValuesTLAD();

  /// \brief fills in \p geovals for all observations in the timeframe (\p t1, \p t2],
  //  \p geovals are interpolated trilinearly from \p state at the nearest gridpoints
  void setTrajectory(const StateQG & state, const util::DateTime & t1,
                     const util::DateTime & t2, GomQG & geovals);
  /// \brief fills in \p geovals for all observations in the timeframe (\p t1, \p t2],
  //  as a tangent-linear operator applied to \p inc
  void fillGeoVaLsTL(const IncrementQG & inc, const util::DateTime & t1,
                     const util::DateTime & t2, GomQG & geovals) const;
  /// \brief fills in \p inc as adjoint operator applied to \p geovals for all
  //  observations in the timeframe (\p t1, \p t2]
  void fillGeoVaLsAD(IncrementQG & inc, const util::DateTime & t1,
                     const util::DateTime & t2, const GomQG & geovals) const;

/// Data
 private:
  void print(std::ostream &) const;
  F90getvalues key_;
};
// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_GETVALUESTLAD_H_
