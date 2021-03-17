/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <memory>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "model/GetValuesQG.h"

#include "oops/util/Logger.h"

#include "model/GeometryQG.h"
#include "model/GomQG.h"
#include "model/LocationsQG.h"
#include "model/StateQG.h"

namespace qg {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
GetValuesQG::GetValuesQG(const GeometryQG &, const LocationsQG & locs,
                         const eckit::Configuration & conf)
    : locs_(locs), conf_(conf)
{
  oops::Log::trace() << "GetValuesQG constructor with config " << conf_ << std::endl;
}

// -----------------------------------------------------------------------------
/// Get state values at observation locations
// -----------------------------------------------------------------------------
void GetValuesQG::fillGeoVaLs(const StateQG & state, const util::DateTime & t1,
                              const util::DateTime & t2, GomQG & gom) const
{
  oops::Log::trace() << "GetValuesQG::fillGeoVaLs start" << std::endl;
  // the below call is an example if one wanted a different interpolation type
  const std::string interpType = conf_.getString("interpolation type", "default");

  if (interpType == "default" || (interpType.compare(0, 8, "default_") == 0)) {
    qg_getvalues_interp_f90(locs_, state.fields().toFortran(), t1, t2, gom.toFortran());
  } else {
    std::string err_message("interpolation type option " + interpType + " not supported");
    throw eckit::BadValue(err_message, Here());
  }
  oops::Log::trace() << "GetValuesQG::fillGeoVaLs done" << std::endl;
}

// -----------------------------------------------------------------------------
void GetValuesQG::print(std::ostream & os) const {
  os << "QG GetValues";
}
// -----------------------------------------------------------------------------

}  // namespace qg
