/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "eckit/mpi/Comm.h"
#include "oops/mpi/mpi.h"

namespace util {

// -----------------------------------------------------------------------------

// Print runtime and memory statistics (only for task 0 by default)
void printRunStats(const std::string &, const bool alltasks = false,
                   const eckit::mpi::Comm & comm = oops::mpi::world());

// -----------------------------------------------------------------------------

}  // namespace util
