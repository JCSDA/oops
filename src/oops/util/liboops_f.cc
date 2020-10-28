/*
 * (C) Copyright 2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/util/liboops_f.h"

namespace oops {

// -----------------------------------------------------------------------------

void liboops_initialise_f() {
  LibOOPS::instance().initialise();
}

// -----------------------------------------------------------------------------

void liboops_finalise_f() {
  LibOOPS::instance().finalise();
}

// -----------------------------------------------------------------------------

}  // namespace oops
