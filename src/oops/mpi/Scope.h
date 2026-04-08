/*
 * (C) Crown Copyright 2025-2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string_view>
#include <string>

namespace eckit::mpi {
class Comm;
}

namespace oops {
namespace mpi {

/// \brief Utility object which sets the default eckit communicator,
///        and resets it when it goes out of scope.
///
/// \details Usage example:
///
///          // -> In MPI_COMM_WORLD space
///          {
///            mpi::Scope _mpiScope("sub");
///            // -> Do things on sub communicator.
///          }
///
///          // -> Returns to MPI_COMM_WORLD
///
///
///          Note that this Scope object relies on the structure of the code stack, for itself
///          to operate as a stack...
///
///          TODO(@mo-joshuacolclough 12/03/26): Whilst Atlas has an `atlas::mpi::Scope`,
///            it has a bug (#186) which was recently fixed. There will need to be a new
///            Atlas tag to incorporate the fix.
class Scope {
 public:
  explicit Scope(const std::string_view nextComm);
  explicit Scope(const eckit::mpi::Comm& nextComm);
  ~Scope();

  Scope(const Scope&)            = delete;
  Scope& operator=(const Scope&) = delete;
 private:
  const std::string startComm_;
};

}  // namespace mpi
}  // namespace oops

