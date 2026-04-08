/*
 * (C) Crown Copyright 2025-2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/mpi/Comm.h"

#include "oops/util/Factory.h"

namespace util {

/// \brief Pure virtual interface for CommRedistribution methods.
/// \details A CommRedistribution performs a redistribution of an Atlas field to
///          and from an MPI subcommunicator split group/"color".
///          This has to be done to repartition a field onto a subset of MPI ranks - `toSubComm`.
///          A copy of the field is sent to each color group in the subcommunicator.
///          Contrarily, `toParentComm` redistributes fields from each MPI color back to the
///          "parent" MPI communicator, as a list of fields where each field is the contents
///          from a given MPI color.
///
///          The source and target functionspaces must have the same grid. However,
///          they can be partitioned differenty as they should be redistributed as needed.
class CommRedistribution {
 public:
  /// \brief Redistribute a field on a larger 'parent' MPI communicator to all members on a
  ///        smaller 'sub' MPI communicator.
  virtual void broadcastToSubMembers(const atlas::Field& sourceParent,
                                     atlas::Field& targetSub) const = 0;

  /// \brief Redistribute a field on a smaller 'sub' MPI comm to a larger 'parent' MPI comm.
  ///        Returns an array of atlas::Fields, where each Field contains the redistributed
  ///        field from a given color. E.g, fields[0] - color 0, fields[1] - color 1, etc.
  virtual std::vector<atlas::Field> gatherToParent(const atlas::Field& sourceSub) const = 0;
};

// Define factory.
                                          // Base class.
using CommRedistributionFactory = Factory<CommRedistribution,
                                          // Constructor arguments.
                                          const eckit::mpi::Comm&, const eckit::mpi::Comm&,
                                          const atlas::FunctionSpace&, const atlas::FunctionSpace&>;

}  // namespace util
