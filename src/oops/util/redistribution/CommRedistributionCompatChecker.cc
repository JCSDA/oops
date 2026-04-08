/*
 * (C) Crown Copyright 2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/redistribution/CommRedistributionCompatChecker.h"

namespace util {


CommRedistributionCompatChecker::CommRedistributionCompatChecker(
    const eckit::mpi::Comm& subComm,
    const eckit::mpi::Comm& parentComm,
    const atlas::FunctionSpace& subFSpace,
    const atlas::FunctionSpace& parentFSpace) :
  subCommName_(subComm.name()),
  parentCommName_(parentComm.name()),
  subFSpace_(subFSpace),
  parentFSpace_(parentFSpace),
  subFSpaceOwnedSize_(countOwned(subFSpace_)),
  parentFSpaceOwnedSize_(countOwned(parentFSpace_)),
  subFSpaceLonLatHash_(getLonLatHash(subFSpace_)),
  parentFSpaceLonLatHash_(getLonLatHash(parentFSpace_))
{
  // Basic checks to confirm these communicators are valid.
  ASSERT(subComm.size() < parentComm.size());
  ASSERT(parentCommName_ == parentFSpace_.mpi_comm());
  ASSERT(subCommName_ == subFSpace_.mpi_comm());

  // Confirm that functionspaces are on the same grid.
  ASSERT_MSG(getGridUid(subFSpace_) == getGridUid(parentFSpace_),
    "CommRedistributionCompatChecker: Mismatch in Grid UID between the subcommunicator "
    "functionspace and parent communicator function space. Both need to be on the same grid.");
}

// ------------------------------------------------------------------------------------------------
namespace {

void compareOwnedGlobalIndices(const atlas::FunctionSpace& fs1, const atlas::FunctionSpace& fs2) {
  // Special case for PointCloud - if both functionspaces are PointCloud skip the check assuming
  // lonlats are checked.
  if (atlas::functionspace::PointCloud(fs1) && atlas::functionspace::PointCloud(fs2)) {
    return;
  } else if (atlas::functionspace::PointCloud(fs1) ^ atlas::functionspace::PointCloud(fs2)) {
    // XOR - If only one of the functionspaces is a PointCloud then they do not match.
    ABORT("compareOwnedGlobalIndices: Only one functionspace is a PointCloud - "
          "functionspace mismatch.");
  }

  // If neither are a PointCloud, continue...
  const auto gidx1 = atlas::array::make_view<atlas::gidx_t, 1>(fs1.global_index());
  const auto ghost1 = atlas::array::make_view<int, 1>(fs1.ghost());
  const auto gidx2 = atlas::array::make_view<atlas::gidx_t, 1>(fs2.global_index());
  const auto ghost2 = atlas::array::make_view<int, 1>(fs2.ghost());

  for (atlas::idx_t ij = 0; ij < gidx1.shape(0); ++ij) {
    // Ensure all ghost points are accounted for. Otherwise mismatch could be missed.
    if (ghost1(ij) == 0 || ghost2(ij) == 0) {
      ASSERT(ghost1(ij) == 0);
      ASSERT(ghost2(ij) == 0);

      ASSERT(gidx1(ij) == gidx2(ij));
    }
  }
}

}  // namespace

// ------------------------------------------------------------------------------------------------

void CommRedistributionCompatChecker::assertFSpaceMatchesSubSpace(
    const atlas::FunctionSpace& fspace) const {
  ASSERT(fspace.mpi_comm() == subCommName_);

  // Note that it doesn't matter if the halos are different sizes, so long as
  // the owned points are the same. Redistributions only concern owned data.
  ASSERT(countOwned(fspace) == subFSpaceOwnedSize_);

  ASSERT_MSG(getLonLatHash(fspace) == subFSpaceLonLatHash_,
             "CommRedistributionCompatChecker::assertFSpaceMatchesSubSpace: Mismatch in "
             "lonlats detected.");

  compareOwnedGlobalIndices(subFSpace_, fspace);
}

// ------------------------------------------------------------------------------------------------

void CommRedistributionCompatChecker::assertFSpaceMatchesParentSpace(
    const atlas::FunctionSpace& fspace) const {
  ASSERT(fspace.mpi_comm() == parentCommName_);
  ASSERT(countOwned(fspace) == parentFSpaceOwnedSize_);

  ASSERT_MSG(getLonLatHash(fspace) == parentFSpaceLonLatHash_,
             "CommRedistributionCompatChecker::assertFSpaceMatchesParentSpace: Mismatch in "
             "lonlats detected.");
  compareOwnedGlobalIndices(parentFSpace_, fspace);
}

// ------------------------------------------------------------------------------------------------


}  // namespace util

