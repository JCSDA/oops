/*
 * (C) Crown Copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/FieldSetSubCommunicators.h"

#include <algorithm>
#include <sstream>
#include <tuple>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/field/for_each.h"

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/Logger.h"

namespace {

// -----------------------------------------------------------------------------

/// Details: The root of MPI communicator `subComm` has rank `subRoot` on this communicator.
///          This function determines the rank of this processing element on the larger
///          communicator `comm`, and a vector of such ranks for all different instances
///          of `subComm`.
std::tuple<std::vector<size_t>, size_t> getLocalRootRanks(
                                            const size_t subRoot,
                                            const eckit::mpi::Comm & subComm,
                                            const eckit::mpi::Comm & comm) {
  // Check local subroot exists
  if (subRoot > subComm.size()) {
    std::stringstream errorMsg;
    errorMsg << "Asked for root " << subRoot << " on sub-communicator "
             << subComm.name() << " of size " << subComm.size() << "."
             << std::endl;
    throw eckit::UserError(errorMsg.str(), Here());
  }

  // Rank of local root on global communicator
  size_t root(0);
  if (subComm.rank() == subRoot) {
    root = comm.rank();
  }

  // Share global root value within sub-communicator
  subComm.broadcast(root, subRoot);

  // Gather local root ranks from all sub-communicators
  std::vector<size_t> roots(comm.size());
  comm.allGather(root, roots.begin(), roots.end());

  // Only keep one instance of each local root
  std::sort(roots.begin(), roots.end());
  auto last = std::unique(roots.begin(), roots.end());
  roots.erase(last, roots.end());

  // Get index of current sub-communicator in vector of `roots`
  const size_t subCommIndex = std::distance(roots.begin(),
                                            std::find(roots.begin(), roots.end(), root));

  return {roots, subCommIndex};
}
}  // namespace

namespace util {
// -----------------------------------------------------------------------------

void redistributeToSubcommunicator(const atlas::FieldSet & fsetIn,
                                   atlas::FieldSet & fsetOut,
                                   const eckit::mpi::Comm & comm,
                                   const eckit::mpi::Comm & subComm,
                                   const atlas::FunctionSpace & fspaceIn,
                                   const atlas::FunctionSpace & fspaceOut) {
  oops::Log::trace() << "redistributeToSubcommunicator starting, from "
                     << comm.name() << " to " << subComm.name() << std::endl;

  const eckit::mpi::Comm & initialDefaultComm = eckit::mpi::comm();

  // Get ranks of local root PEs on global communicator and sub-communicators
  const size_t subRoot(0);
  auto[roots, subCommIndex] = getLocalRootRanks(subRoot, subComm, comm);

  // Redistribute fields one by one to limit memory usage
  for (const auto & fieldIn : fsetIn) {
    const auto fieldName = fieldIn.name();
    eckit::mpi::setCommDefault(comm.name().c_str());

    // TODO(Mayeul): add different data types (int, long, float)
    ASSERT(fieldIn.datatype() == atlas::array::DataType::kind<double>());

    // 1. Gather on each sub-communicator
    std::vector<atlas::Field> rootFields;
    for (const size_t otherRoot : roots) {
      rootFields.push_back(fspaceIn.createField<double>(
                  atlas::option::name(fieldName)
                | atlas::option::levels(fieldIn.levels())
                | atlas::option::global(otherRoot)));

      fspaceIn.gather(fieldIn, rootFields.back());
    }

    // 2. Copy global fields defined on global communicator to global field
    //    defined on sub-communicator
    eckit::mpi::setCommDefault(subComm.name().c_str());
    // Changing the default communicator means we should change the owner of the global field
    // from root to sub-root. Unfortunately, atlas field offers no way to change this.
    // As an alternative, we copy the field into another global field on the same PE
    // with the correct root numbering.
    atlas::Field rootField = fspaceOut.createField<double>(
                atlas::option::name(fieldName)
              | atlas::option::levels(fieldIn.levels())
              | atlas::option::global(subRoot));
    atlas::array::make_view<double, 2>(rootField).assign(
                atlas::array::make_view<double, 2>(rootFields[subCommIndex]));

    // 3. Scatter from local root PEs to each subcommunicator
    if (!fsetOut.has(fieldName)) {
      fsetOut.add(fspaceOut.createField<double>(
                         atlas::option::name(fieldName)
                       | atlas::option::levels(fieldIn.levels())));
    }

    fspaceOut.scatter(rootField, fsetOut[fieldName]);

    fsetOut[fieldName].set_dirty(fieldIn.dirty());
  }

  eckit::mpi::setCommDefault(initialDefaultComm.name().c_str());
  oops::Log::trace() << "redistributeToSubcommunicator done, from "
                     << comm.name() << " to " << subComm.name()
                     << std::endl;
}

// -----------------------------------------------------------------------------

void gatherAndSumFromSubcommunicator(const atlas::FieldSet & fsetIn,
                                     atlas::FieldSet & fsetOut,
                                     const eckit::mpi::Comm & subComm,
                                     const eckit::mpi::Comm & comm,
                                     const atlas::FunctionSpace & fspaceIn,
                                     const atlas::FunctionSpace & fspaceOut) {
  oops::Log::trace() << "gatherAndSumFromSubcommunicator starting, from "
                     << subComm.name() << " to " << comm.name()
                     << std::endl;
  const eckit::mpi::Comm & initialDefaultComm = eckit::mpi::comm();

  const size_t subRoot(0);
  auto[roots, subCommIndex] = getLocalRootRanks(subRoot, subComm, comm);

  // Redistribute fields one by one to limit memory usage
  for (const auto & fieldIn : fsetIn) {
    const auto fieldName = fieldIn.name();

    // 1. Gather field on root PE of subcommunicator
    eckit::mpi::setCommDefault(subComm.name().c_str());

    // TODO(Mayeul): add different data types (int, long, float)
    ASSERT(fieldIn.datatype() == atlas::array::DataType::kind<double>());

    atlas::Field rootField = fspaceIn.createField<double>(
                  atlas::option::name(fieldName)
                | atlas::option::levels(fieldIn.levels())
                | atlas::option::global(subRoot));
    fspaceIn.gather(fieldIn, rootField);

    // 2. Copy sub-communicator global field to global-communicator global field
    eckit::mpi::setCommDefault(comm.name().c_str());
    std::vector<atlas::Field> rootFields;
    for (const size_t otherRoot : roots) {
      rootFields.push_back(fspaceOut.createField<double>(
                  atlas::option::name(fieldName)
                | atlas::option::levels(fieldIn.levels())
                | atlas::option::global(otherRoot)));
    }

    atlas::array::make_view<double, 2>(rootFields[subCommIndex]).assign(
                atlas::array::make_view<double, 2>(rootField));

    // 3. Scatter global field from root PE to whole communicator and sum

    // Initialize field with zeros
    if (!fsetOut.has(fieldName)) {
      fsetOut.add(fspaceOut.createField<double>(
                               atlas::option::name(fieldName)
                             | atlas::option::levels(fieldIn.levels())));
    }
    auto & fieldOut = fsetOut[fieldName];
    atlas::array::make_view<double, 2>(fieldOut).assign(0.0);

    auto fieldOutTmp = fspaceOut.createField<double>(atlas::option::name(fieldName)
                                          | atlas::option::levels(fieldOut.levels()));
    for (const auto & globalRootField : rootFields) {
      // Scatter
      fspaceOut.scatter(globalRootField, fieldOutTmp);

      // Sum
      atlas::field::for_each_value(atlas::execution::par_unseq,
                                   fieldOutTmp, fieldOut,
                                   [&](const double a, double& b) { b += a; });
    }

    // 4. Set dirty halos if any of the input fields is dirty
    int dirty = static_cast<int>(fieldIn.dirty());
    comm.allReduceInPlace(dirty, eckit::mpi::max());
    fieldOut.set_dirty(static_cast<bool>(dirty));
  }

  eckit::mpi::setCommDefault(initialDefaultComm.name().c_str());
  oops::Log::trace() << "gatherAndSumFromSubcommunicator done, from "
                     << subComm.name() << " to " << comm.name()
                     << std::endl;
}

// -----------------------------------------------------------------------------
}  // namespace util
