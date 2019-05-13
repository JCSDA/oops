/*
 * (C) Copyright 2017-2018 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "model/qg_obsvec_random_f.h"
#include "model/ObsSpaceQG.h"

// -----------------------------------------------------------------------------
namespace qg {
// -----------------------------------------------------------------------------

void qg_obsvec_random_f(const ObsSpaceQG * odb, const int & nn, double * xx) {
  odb->random(nn, xx);
}

// -----------------------------------------------------------------------------

}  // namespace qg
