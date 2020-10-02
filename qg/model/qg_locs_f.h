/*
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QG_MODEL_QG_LOCS_F_H_
#define QG_MODEL_QG_LOCS_F_H_

#include "atlas/field.h"
#include "atlas/functionspace/PointCloud.h"

#include "oops/qg/LocationsQG.h"
#include "oops/util/DateTime.h"

// ------------------------------------------------------------------------------
// These functions provide tools for interfacing Fortran and C++ string objects
// ------------------------------------------------------------------------------

namespace qg {

extern "C" {
  int qg_locs_nlocs_f90(qg::LocationsQG*);
  atlas::field::FieldImpl* qg_locs_lonlat_f90(qg::LocationsQG*);
  atlas::field::FieldImpl* qg_locs_altitude_f90(qg::LocationsQG*);
  util::DateTime& qg_locs_times_f90(qg::LocationsQG*, size_t &);
}

}  // namespace qg

#endif  // QG_MODEL_QG_LOCS_F_H_
