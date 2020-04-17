/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QG_MODEL_GETVALUESQG_H_
#define QG_MODEL_GETVALUESQG_H_

#include <ostream>
#include <string>

#include "model/QgFortran.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace qg {
  class GomQG;
  class LocationsQG;
  class GeometryQG;
  class StateQG;

/// \brief used for getting state values at observation locations
// -----------------------------------------------------------------------------
class GetValuesQG : public util::Printable,
                    private util::ObjectCounter<GetValuesQG> {
 public:
  static const std::string classname() {return "qg::GetValuesQG";}

/// \brief saves all locations \p locs to use during filling GeoVaLs
  GetValuesQG(const GeometryQG &, const LocationsQG & locs);
  virtual ~GetValuesQG();

/// \brief fills in \p geovals for all observations in the timeframe (\p t1, \p t2],
/// \p geovals are interpolated trilinearly from \p state at the nearest gridpoints
  void fillGeoVaLs(const StateQG &, const util::DateTime & t1,
                   const util::DateTime & t2, GomQG &) const;

 private:
  void print(std::ostream &) const;
  F90getvalues key_;
};
// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_GETVALUESQG_H_
