/*
 * (C) Crown Copyright 2025-2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "eckit/mpi/Comm.h"

#include "oops/mpi/Scope.h"
#include "oops/util/Logger.h"

namespace oops {
namespace mpi {

Scope::Scope(const std::string_view nextComm) :
  startComm_(eckit::mpi::comm().name())
{
  oops::Log::trace() << "oops::mpi::Scope begun with '" << startComm_
                     << "' { '" << nextComm << "' }" << std::endl;
  eckit::mpi::setCommDefault(nextComm);
}

Scope::Scope(const eckit::mpi::Comm& nextComm) :
  Scope(nextComm.name())
{}

Scope::~Scope() {
  oops::Log::trace() << "oops::mpi::Scope finished, returning to '" << startComm_
                     << "'." << std::endl;
  eckit::mpi::setCommDefault(startComm_);
}

}  // namespace mpi
}  // namespace oops

