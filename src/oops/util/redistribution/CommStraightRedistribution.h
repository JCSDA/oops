/*
 * (C) Crown Copyright 2025-2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "atlas/functionspace.h"

#include "oops/mpi/ColorInfo.h"
#include "oops/util/redistribution/CommRedistribution.h"
#include "oops/util/redistribution/CommRedistributionCompatChecker.h"
#include "oops/util/redistribution/detail/AllToAllRouting.h"

namespace util {

// ------------------------------------------------------------------------------------------------
/// \brief Atlas Field 'straight' redistribution across communicators.
/// \details This redistribution works by creating a direct, global index mapping between ranks
///          on the `parent` communicator, for each color in an MPI split communicator.
///          The mapping is then used to perform an MPI `allToAllv` operation to redistribute
///          an Atlas Field.
///          Once constructed this mapping can be reused for every redistribution.
class CommStraightRedistribution : public CommRedistribution {
 public:
  CommStraightRedistribution(const eckit::mpi::Comm& sub,
                             const eckit::mpi::Comm& parent,
                             const atlas::FunctionSpace& subFSpace,
                             const atlas::FunctionSpace& parentFSpace);

  void broadcastToSubMembers(const atlas::Field& sourceParent,
                             atlas::Field& targetSub) const override;

  std::vector<atlas::Field> gatherToParent(const atlas::Field& sourceSub) const override;

 private:
  static const std::string classname() { return "CommStraightRedistribution"; }

  /// \brief Transfer a field from source to target, on a given MPI color group.
  void transferFieldOnColor(const atlas::Field& source,
                             const size_t color,
                             atlas::Field& target) const;

  const eckit::mpi::Comm& subComm_;
  const eckit::mpi::Comm& parentComm_;
  const oops::mpi::ColorInfo colorInfo_;

  const atlas::FunctionSpace parentFSpace_;
  const CommRedistributionCompatChecker checker_;

  // Mapping from global index -> local index on the local parent partition.
  std::unordered_map<atlas::uidx_t, atlas::idx_t> parentGlobalToLocalIndexMap_;
  // Mapping from global index -> local index on the local subcommunicator partition.
  std::unordered_map<atlas::uidx_t, atlas::idx_t> subGlobalToLocalIndexMap_;

  // color[parent rank[global index list ...]]
  // Lists of global indices to send to each parent rank in the sub -> parent direction.
  std::vector<std::vector<std::vector<atlas::uidx_t>>> routingToParent_;
  // Lists of global indices to recieve on each parent rank in the sub -> parent direction.
  std::vector<std::vector<std::vector<atlas::uidx_t>>> routingToSub_;

  // Utility for an MPI allToAllv operation to hold horizontal send & recv buffer sizes
  // for each rank, and calculate buffer displacements.
  // Must be transformed depending on the number of field levels.
  std::vector<detail::AllToAllRouting> horizontalAllToAllToParent_;

  size_t totalSendCountToParent_, totalRecvCountToParent_;
};

}  // namespace util

