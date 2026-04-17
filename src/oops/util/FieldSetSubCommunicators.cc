/*
 * (C) Crown Copyright 2024-2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/FieldSetSubCommunicators.h"

#include <algorithm>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/field/for_each.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/FieldSetOperations.h"
#include "oops/util/Logger.h"
#include "oops/util/redistribution/CommStraightRedistribution.h"

namespace util {

// -----------------------------------------------------------------------------

void redistributeToSubcommunicator(const CommRedistribution & redist,
                                   const atlas::FieldSet & fsetIn,
                                   atlas::FieldSet & fsetOut,
                                   const atlas::FunctionSpace & fspaceOut) {
  if (fsetIn.size() <= 0) return;

  // Ensure that field names are sorted before performing redistribution.
  // Otherwise there could be mix-ups where colors have differing orderings.
  std::vector<std::string> sortedFieldNames = fsetIn.field_names();
  std::sort(sortedFieldNames.begin(), sortedFieldNames.end());

  // Redistribute fields one by one to limit memory usage
  for (const std::string& fieldName : sortedFieldNames) {
    const atlas::Field& fieldIn = fsetIn[fieldName];

    if (!fsetOut.has(fieldName)) {
      fsetOut.add(fspaceOut.createField<double>(
                         atlas::option::name(fieldName)
                       | atlas::option::levels(fieldIn.levels())));
    }

    redist.broadcastToSubMembers(fieldIn, fsetOut[fieldName]);
  }
}

// -----------------------------------------------------------------------------

void redistributeToSubcommunicator(const atlas::FieldSet & fsetIn,
                                   atlas::FieldSet & fsetOut,
                                   const eckit::mpi::Comm & comm,
                                   const eckit::mpi::Comm & subComm,
                                   const atlas::FunctionSpace & fspaceIn,
                                   const atlas::FunctionSpace & fspaceOut) {
  // Create a straight redistribution
  const CommStraightRedistribution redist(subComm, comm, fspaceOut, fspaceIn);
  redistributeToSubcommunicator(redist, fsetIn, fsetOut, fspaceOut);
}

// -----------------------------------------------------------------------------

void gatherAndSumFromSubcommunicator(const CommRedistribution & redist,
                                     const atlas::FieldSet & fsetIn,
                                     atlas::FieldSet & fsetOut,
                                     const eckit::mpi::Comm & comm,
                                     const atlas::FunctionSpace & fspaceIn,
                                     const atlas::FunctionSpace & fspaceOut) {
  // Ensure that field names are sorted before performing redistribution.
  // Otherwise there could be mix-ups where colors have differing orderings.
  std::vector<std::string> sortedFieldNames = fsetIn.field_names();
  std::sort(sortedFieldNames.begin(), sortedFieldNames.end());

  // Redistribute fields one by one to limit memory usage
  for (const std::string& fieldName : sortedFieldNames) {
    const atlas::Field& fieldIn = fsetIn[fieldName];

    // Initialize field with zeros
    if (!fsetOut.has(fieldName)) {
      fsetOut.add(fspaceOut.createField<double>(
                               atlas::option::name(fieldName)
                             | atlas::option::levels(fieldIn.levels())));
    }
    auto & fieldOut = fsetOut[fieldName];
    atlas::array::make_view<double, 2>(fieldOut).assign(0.0);

    const std::vector<atlas::Field> subCommFields = redist.gatherToParent(fieldIn);

    for (const auto & subCommField : subCommFields) {
      ASSERT(subCommField.functionspace().mpi_comm() == fspaceOut.mpi_comm());

      // Sum
      atlas::field::for_each_value(atlas::execution::par_unseq,
                                   subCommField, fieldOut,
                                   [&](const double a, double& b) { b += a; });
    }

    // 4. Set dirty halos if any of the input fields is dirty
    fieldOut.set_dirty();
  }
}

// -----------------------------------------------------------------------------

void gatherAndSumFromSubcommunicator(const atlas::FieldSet & fsetIn,
                                     atlas::FieldSet & fsetOut,
                                     const eckit::mpi::Comm & subComm,
                                     const eckit::mpi::Comm & comm,
                                     const atlas::FunctionSpace & fspaceIn,
                                     const atlas::FunctionSpace & fspaceOut) {
  // Create a straight redistribution
  const CommStraightRedistribution redist(subComm, comm, fspaceIn, fspaceOut);
  gatherAndSumFromSubcommunicator(redist, fsetIn, fsetOut, comm, fspaceIn, fspaceOut);
}

// -----------------------------------------------------------------------------
}  // namespace util
