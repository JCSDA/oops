/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/workflow.h"

#include <cstdlib>

#include "eckit/exception/Exceptions.h"
#include "oops/mpi/mpi.h"
#include "oops/util/Logger.h"

namespace util {

// -----------------------------------------------------------------------------

static bool has_ecflow = false;

// -----------------------------------------------------------------------------

void use_ecflow() {has_ecflow = true;}

// -----------------------------------------------------------------------------

void update_workflow_meter(const std::string & name, const int value) {
  if (has_ecflow) {
    if (oops::mpi::world().rank() == 0) {
      const std::string ecfsignal = "ecflow_client --meter=" + name + " " + std::to_string(value);
      const char * cmd = ecfsignal.c_str();
      const int ret = std::system(cmd);
      ASSERT(ret == 0);
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace util
