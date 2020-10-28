/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "atlas/field.h"
#include "atlas/functionspace/PointCloud.h"

#include "model/LocationsQG.h"
#include "model/qg_locs_f.h"
#include "oops/util/DateTime.h"

namespace qg {

// -----------------------------------------------------------------------------
int qg_locs_nlocs_f90(qg::LocationsQG* locs) {
    return locs->size();
}
atlas::field::FieldImpl* qg_locs_lonlat_f90(qg::LocationsQG* locs) {
    return locs->lonlat().get();
}
atlas::field::FieldImpl* qg_locs_altitude_f90(qg::LocationsQG* locs) {
    return locs->altitude().get();
}
util::DateTime& qg_locs_times_f90(qg::LocationsQG* locs, size_t & idx) {
    return locs->times(idx);
}
// -----------------------------------------------------------------------------

}  // namespace qg
