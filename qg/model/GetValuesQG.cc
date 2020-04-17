/*
 * (C) Copyright 2019-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

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
GetValuesQG::GetValuesQG(const GeometryQG & geom, const LocationsQG & locs) {
  qg_getvalues_create_f90(key_, geom.toFortran(), locs.toFortran());
}
// -----------------------------------------------------------------------------
GetValuesQG::~GetValuesQG() {
  qg_getvalues_delete_f90(key_);
}
// -----------------------------------------------------------------------------
/// Get state values at observation locations
// -----------------------------------------------------------------------------
void GetValuesQG::fillGeoVaLs(const StateQG & state, const util::DateTime & t1,
                              const util::DateTime & t2, GomQG & gom) const {
  qg_getvalues_interp_f90(key_, state.fields().toFortran(),
                          t1, t2, gom.toFortran());
}
// -----------------------------------------------------------------------------
void GetValuesQG::print(std::ostream & os) const {
  os << "GetValues" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace qg
