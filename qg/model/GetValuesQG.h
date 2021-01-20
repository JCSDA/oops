/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QG_MODEL_GETVALUESQG_H_
#define QG_MODEL_GETVALUESQG_H_

#include <memory>
#include <ostream>
#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "oops/qg/LocationsQG.h"
#include "oops/qg/QgFortran.h"

namespace eckit {
  class Configuration;
}

namespace qg {
  class GomQG;
  class GeometryQG;
  class StateQG;

/// \brief used for getting state values at observation locations
// -----------------------------------------------------------------------------
class GetValuesQG : public util::Printable,
                    private util::ObjectCounter<GetValuesQG> {
 public:
  static const std::string classname() {return "qg::GetValuesQG";}

/// \brief saves all locations \p locs to use during filling GeoVaLs
  GetValuesQG(const GeometryQG &, const LocationsQG & locs,
              const eckit::Configuration &);
  ~GetValuesQG() {}

/// \brief fills in \p geovals for all observations in the timeframe (\p t1, \p t2],
/// \p geovals are interpolated trilinearly from \p state at the nearest gridpoints
  void fillGeoVaLs(const StateQG &, const util::DateTime & t1,
                   const util::DateTime & t2, GomQG &) const;

 private:
  void print(std::ostream &) const;
  LocationsQG locs_;
  eckit::LocalConfiguration conf_;
};
// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_GETVALUESQG_H_
