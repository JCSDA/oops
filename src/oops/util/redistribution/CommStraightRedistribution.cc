/*
 * (C) Crown Copyright 2025-2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <optional>
#include <algorithm>

#include "eckit/exception/Exceptions.h"

#include "oops/mpi/mpi.h"
#include "oops/mpi/Scope.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/defines.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/redistribution/CommStraightRedistribution.h"
#include "oops/util/Timer.h"

namespace util {

// ------------------------------------------------------------------------------------------------

// Register in factory.
static CommRedistributionFactory::Maker<CommStraightRedistribution>
    makerStraightRedist_("straight");

// ------------------------------------------------------------------------------------------------
namespace {

// ------------------------------------------------------------------------------------------------

// This function creates a mapping for the whole grid of what points are owned by what partition,
// in an MPI subcommunicator, for a given color.
std::unordered_map<atlas::uidx_t, size_t> getGlobalIndexOwners(
    const eckit::mpi::Comm& parentComm,
    const atlas::FunctionSpace& subFSpace,
    const std::optional<bool> isActiveColor = std::nullopt) {
  // In the case that the subject is a functionspace on a sub communicator,
  // each color must be treated differently. Only if this is the active color then
  // global indices are added. Otherwise, this rank does not own any of the field!
  // Each MPI color group will own completely separate field data.
  //
  // i.e If the subject is not a subComm, then add these global indices.
  //     OR, if the subject is a subComm, then ONLY when on the active color add the
  //     global indices. (isActiveColor contains has a value -> is a subcomm, and it's true).
  //
  bool shouldAddIndicesOnThisRank = true;
  if (isActiveColor.has_value()) {
    shouldAddIndicesOnThisRank = *isActiveColor;
  }

  // Get number of points expected to be received, across parent communicator!
  // Will equal 0 for those partitions not in this color group.
  std::vector<int> partitionSizes(parentComm.size(), -1);

  const size_t numOwned = shouldAddIndicesOnThisRank ? countOwned(subFSpace) : 0;
  parentComm.allGather(static_cast<int>(numOwned), partitionSizes.begin(), partitionSizes.end());

  const size_t totalNumPoints = std::accumulate(partitionSizes.begin(), partitionSizes.end(), 0);

  // Gather global indices.
  std::vector<atlas::uidx_t> globalIndices;
  globalIndices.reserve(totalNumPoints);

  if (shouldAddIndicesOnThisRank) {
    const auto gidx = atlas::array::make_view<const atlas::uidx_t, 1>(subFSpace.global_index());
    const auto ghost = atlas::array::make_view<const int, 1>(subFSpace.ghost());

    for (atlas::idx_t ij = 0; ij < gidx.shape(0); ++ij) {
      if (ghost(ij) == 0) {
        globalIndices.emplace_back(gidx(ij));
      }
    }
  }

  ASSERT(globalIndices.size() == numOwned);
  // Gather
  oops::mpi::allGatherv(parentComm, globalIndices);
  ASSERT(globalIndices.size() == totalNumPoints);

  std::unordered_map<atlas::uidx_t, size_t> gidxOwners;
  gidxOwners.reserve(totalNumPoints);

  auto buf = globalIndices.begin();
  for (size_t rank = 0; rank < parentComm.size(); ++rank) {
    const size_t rankSize = partitionSizes[rank];
    for (size_t i = 0; i < rankSize; ++i) {
      auto[_it, success] = gidxOwners.emplace(*buf++, rank);
      ASSERT(success);  // Already existed in map... Overwritten!
    }
  }

  ASSERT(gidxOwners.size() == totalNumPoints);

  return gidxOwners;
}

}  // namespace

// ------------------------------------------------------------------------------------------------

CommStraightRedistribution::CommStraightRedistribution(const eckit::mpi::Comm& sub,
                                                       const eckit::mpi::Comm& parent,
                                                       const atlas::FunctionSpace& subFSpace,
                                                       const atlas::FunctionSpace& parentFSpace) :
  subComm_(sub),
  parentComm_(parent),
  colorInfo_(subComm_, parentComm_),
  parentFSpace_(parentFSpace),
  checker_(sub, parent, subFSpace, parentFSpace)
{
  Timer _timer(classname(), "CommStraightRedistribution");
  oops::Log::trace() << "CommStraightRedistribution::CommStraightRedistribution start."
                     << std::endl;

  // Build parent ownership map.
  const std::unordered_map<atlas::uidx_t, size_t> parentOwnershipMap =
      getGlobalIndexOwners(parentComm_, parentFSpace_);

  // Build color ownership maps.
  std::vector<std::unordered_map<atlas::uidx_t, size_t>> colorOwnershipMaps;
  colorOwnershipMaps.reserve(colorInfo_.numColors());

  for (size_t col = 0; col < colorInfo_.numColors(); ++col) {
    const bool isActiveColor = col == colorInfo_.color();

    colorOwnershipMaps.emplace_back(
        getGlobalIndexOwners(parentComm_,
                             subFSpace,
                             isActiveColor));

    // Confirm that the color mappings contain the same global indices
    // as the parent mappings.
    if constexpr (oops::OOPS_BUILD_TYPE_DEBUG) {
      const auto& subMap = colorOwnershipMaps.back();

      if (subMap.size() != parentOwnershipMap.size()) {
        ABORT("Sub functionspace global index ownership mapping is not the same size as the parent "
              "functionspace global index ownership mapping. parent {" +
              std::to_string(parentOwnershipMap.size()) +
              "} != {" + std::to_string(subMap.size()) +
              "} sub");
      }

      bool passed = true;
      for (const auto& [gidxParent, _parentOwner] : parentOwnershipMap) {
        if (subMap.find(gidxParent) == subMap.end()) {
          oops::Log::error() << "Sub map missing gidx => " << gidxParent << std::endl;
          passed = false;
        }
      }
      if (!passed) {
        ABORT("Inconsistency between sub functionspace global ownership map and the parent "
              "global ownership map.");
      }
    }
  }

  // Build index mappings from global -> local.
  const auto parentGidx =
      atlas::array::make_view<const atlas::uidx_t, 1>(parentFSpace.global_index());
  const auto parentGhost = atlas::array::make_view<int, 1>(parentFSpace.ghost());
  parentGlobalToLocalIndexMap_.reserve(parentGidx.shape(0));
  for (atlas::idx_t ij = 0; ij < parentGidx.shape(0); ++ij) {
    if (parentGhost(ij) == 0) {
      auto[_it, success] = parentGlobalToLocalIndexMap_.emplace(parentGidx(ij), ij);
      ASSERT(success);   // Already in map... Overwritten!
    }
  }

  const auto subGidx =
      atlas::array::make_view<const atlas::uidx_t, 1>(subFSpace.global_index());
  subGlobalToLocalIndexMap_.reserve(subGidx.shape(0));
  const auto subGhost = atlas::array::make_view<int, 1>(subFSpace.ghost());
  for (atlas::idx_t ij = 0; ij < subGidx.shape(0); ++ij) {
    if (subGhost(ij) == 0) {
      auto[_it, success] = subGlobalToLocalIndexMap_.emplace(subGidx(ij), ij);
      ASSERT(success);   // Already in map... Overwritten!
    }
  }

  // Build send/recv vectors.
  const size_t ownedParentSize = countOwned(parentFSpace);
  const size_t ownedSubSize = countOwned(subFSpace);
  const size_t numColors = colorInfo_.numColors();

  routingToParent_.reserve(numColors);
  routingToSub_.reserve(numColors);

  for (size_t color = 0; color < numColors; ++color) {
    // Global indices to send to each rank.
    std::vector<std::vector<atlas::uidx_t>> sendsOnColor;
    std::vector<std::vector<atlas::uidx_t>> recvsOnColor;

    for (size_t rank = 0; rank < parentComm_.size(); ++rank) {
      sendsOnColor.emplace_back();
      recvsOnColor.emplace_back();
      // Possible optimisation: overestimate vector sizes as size of this PE.
    }
    std::vector<int> sendCountsOnColor(parentComm_.size(), 0);
    std::vector<int> recvCountsOnColor(parentComm_.size(), 0);

    // For each global index (& across whole grid), for this particular color group
    // figure out what needs to be sent to where, and received from where,
    // for the direction going sub -> parent.
    for (const auto& [gidx, parentOwner] : parentOwnershipMap) {
      const size_t subOwner = colorOwnershipMaps.at(color).at(gidx);

      // Send
      if (subOwner == parentComm_.rank()) {
        sendsOnColor[parentOwner].emplace_back(gidx);
        ++sendCountsOnColor[parentOwner];
      }
      // Recv
      if (parentOwner == parentComm_.rank()) {
        recvsOnColor[subOwner].emplace_back(gidx);
        ++recvCountsOnColor[subOwner];
      }
    }

    routingToParent_.push_back(sendsOnColor);
    routingToSub_.push_back(recvsOnColor);

    totalRecvCountToParent_ = std::accumulate(recvCountsOnColor.begin(),
                                              recvCountsOnColor.end(),
                                              0);
    ASSERT(totalRecvCountToParent_ == ownedParentSize);

    totalSendCountToParent_ = std::accumulate(sendCountsOnColor.begin(),
                                              sendCountsOnColor.end(),
                                              0);
    const bool isActiveColor = colorInfo_.color() == color;
    ASSERT(totalSendCountToParent_ == (isActiveColor ? ownedSubSize : 0));

    horizontalAllToAllToParent_.emplace_back(sendCountsOnColor, recvCountsOnColor);
  }


  oops::Log::trace() << "CommStraightRedistribution::CommStraightRedistribution finished."
                     << std::endl;
}

// ------------------------------------------------------------------------------------------------

void CommStraightRedistribution::transferFieldOnColor(const atlas::Field& source,
                                                      const size_t color,
                                                      atlas::Field& target) const {
  oops::Log::trace() << "CommStraightRedistribution::transferFieldOnColor starting" << std::endl;

  ASSERT(source.levels() == target.levels());

  enum Direction { eToParent, eToSub };
  Direction direction = eToParent;
  if (source.functionspace().mpi_comm() == parentComm_.name()) {
    // Confirm that each field's functionspace matches the ones
    // this redistribution is setup for.
    checker_.assertFieldMatchesSubSpace(target);
    checker_.assertFieldMatchesParentSpace(source);
    direction = eToSub;
  } else if (source.functionspace().mpi_comm() == subComm_.name()) {
    checker_.assertFieldMatchesSubSpace(source);
    checker_.assertFieldMatchesParentSpace(target);
    direction = eToParent;
  } else {
    ABORT("Fields in transferFieldOnColor should be on different MPI communicators.");
  }

  // Write the source field to a std::vector buffer in order of rank.
  const auto sourceView = atlas::array::make_view<const double, 2>(source);

  std::vector<double> sendBuffer;
  const size_t sendBufferHorizontalSize = countOwned(source.functionspace());
  // Confirm send field has same size.
  ASSERT(countOwned(source.functionspace()) == sendBufferHorizontalSize);
  sendBuffer.reserve(sendBufferHorizontalSize * source.levels());

  // NOTE: Direction switches what is being sent and what is being receieved.
  //       The mappings stored correspond to what would be sent/received in the gatherToParent
  //       direction. In the case of `broadcastToSubMembers`, these need switching.
  const auto& localIndexMappingSend = direction == eToParent ? subGlobalToLocalIndexMap_ :
                                                               parentGlobalToLocalIndexMap_;
  // Confirm global indices match
  if constexpr (oops::OOPS_BUILD_TYPE_DEBUG) {
    const auto gidxSource = atlas::array::make_view<atlas::uidx_t, 1>(
        source.functionspace().global_index());

    size_t checked = 0;
    for (const auto& [gidx, localidx] : localIndexMappingSend) {
      if (gidxSource(localidx) != gidx) {
        oops::Log::error() << "CommStraightRedistribution::transferFieldOnColor: "
                              "Global index in source field does not match for local index '"
                           << localidx << "'... " << gidx << " != " << gidxSource(localidx)
                           << ". Checked " << checked << " indices." << std::endl;
        ABORT("CommStraightRedistribution::transferFieldOnColor failed: Source field gidx"
              " mismatch.");
      }
      ++checked;
    }
  }

  const auto& sends = direction == eToParent ? routingToParent_[color] :
                                               routingToSub_[color];
  for (const auto& sendToRank : sends) {
    for (const atlas::uidx_t gidx : sendToRank) {
      // Find matching gidx in local index
      const atlas::idx_t localIdx = localIndexMappingSend.at(gidx);

      for (atlas::idx_t k = 0; k < sourceView.shape(1); ++k) {
        sendBuffer.emplace_back(sourceView(localIdx, k));
      }
    }
  }
  // ------

  // Setup recieve buffer.
  const auto& localIndexMappingRecv = direction == eToParent ? parentGlobalToLocalIndexMap_ :
                                                               subGlobalToLocalIndexMap_;
  const size_t recvOwnedSize = localIndexMappingRecv.size();
  ASSERT(countOwned(target.functionspace()) == recvOwnedSize);

  if constexpr (oops::OOPS_BUILD_TYPE_DEBUG) {
    const auto gidxTarget = atlas::array::make_view<atlas::uidx_t, 1>(
        target.functionspace().global_index());

    size_t checked = 0;
    for (const auto& [gidx, localidx] : localIndexMappingRecv) {
      if (gidxTarget(localidx) != gidx) {
        oops::Log::error() << "CommStraightRedistribution::transferFieldOnColor: "
                              "Global index in target field does not match for local index '"
                           << localidx << "'... " << gidx << " != " << gidxTarget(localidx)
                           << ". Checked " << checked << " indices." << std::endl;
        ABORT("CommStraightRedistribution::transferFieldOnColor failed: Target field gidx"
              " mismatch.");
      }
      ++checked;
    }
  }

  const size_t fieldLevels = target.levels();

  constexpr double kFillValue = -9e9;
  std::vector<double> recvBuffer(recvOwnedSize * fieldLevels, kFillValue);

  // Scale horizontal routing displacements by number of levels.
  const detail::AllToAllRouting routing = [&]() {
    switch (direction) {
      case eToParent:
        return horizontalAllToAllToParent_[color] * fieldLevels;
      case eToSub:
        return horizontalAllToAllToParent_[color].invert() * fieldLevels;
    }
    // Unreachable. Surpresses compiler warning.
    return detail::AllToAllRouting();
  }();

  // ALL TO ALL
  routing.execute(parentComm_, sendBuffer, recvBuffer);

  // Write received buffer to the target Field in expected recieve order.
  const auto& recvs = direction == eToParent ? routingToSub_[color] :
                                               routingToParent_[color];
  auto outFieldView = atlas::array::make_view<double, 2>(target);

  auto bufferIter = recvBuffer.begin();
  for (const std::vector<atlas::uidx_t>& recvsFromRank : recvs) {
    for (const atlas::uidx_t gidxRecieved : recvsFromRank) {
      const atlas::idx_t localIdx = localIndexMappingRecv.at(gidxRecieved);

      // Read from buffer
      for (atlas::idx_t k = 0; k < outFieldView.shape(1); ++k) {
        if constexpr (oops::OOPS_BUILD_TYPE_DEBUG) {   // Check only on debug
          if (*bufferIter == kFillValue) {
            const std::string err = "Global index " + std::to_string(gidxRecieved) +
                                    " not populated on rank " +
                                    std::to_string(parentComm_.rank()) + ".";
            ABORT(err);
          }
        }

        outFieldView(localIdx, k) = *bufferIter++;
      }
    }
  }

  target.set_dirty();
  oops::Log::trace() << "CommStraightRedistribution::transferFieldOnColor finished" << std::endl;
}

// ------------------------------------------------------------------------------------------------

std::vector<atlas::Field> CommStraightRedistribution::gatherToParent(
    const atlas::Field& sourceSub) const {
  Timer _timer(classname(), "gatherToParent");

  oops::Log::trace() << "CommStraightRedistribution::gatherToParent starting, from "
                     << subComm_.name() << " to " << parentComm_.name() << std::endl;

  const oops::mpi::Scope _mpiScope(parentComm_);

  ASSERT_MSG(sourceSub.functionspace().mpi_comm() == subComm_.name(),
             "CommStraightRedistribution::gatherToParent: Subcomm field does not "
             "have the same MPI communicator as the sub MPI communicator this redistribution has "
             "been set up with. Has the field been set up correctly?");

  const auto fieldConf = atlas::option::name(sourceSub.name()) |
                         atlas::option::levels(sourceSub.levels()) |
                         atlas::option::datatype(sourceSub.datatype());

  const size_t numColors = colorInfo_.numColors();

  std::vector<atlas::Field> parentFields;

  for (size_t color = 0; color < numColors; ++color) {
    atlas::Field outField = parentFSpace_.createField(fieldConf);

    transferFieldOnColor(sourceSub, color, outField);

    parentFields.push_back(outField);
  }

  oops::Log::trace() << "CommStraightRedistribution::gatherToParent finished, from "
                     << subComm_.name() << " to " << parentComm_.name() << std::endl;
  return parentFields;
}

// ------------------------------------------------------------------------------------------------

void CommStraightRedistribution::broadcastToSubMembers(const atlas::Field& sourceParent,
                                                       atlas::Field& targetSub) const {
  Timer _timer(classname(), "broadcastToSubMembers");
  oops::Log::trace() << "CommStraightRedistribution::broadcastToSubMembers starting, from "
                     << parentComm_.name() << " to " << subComm_.name() << std::endl;

  const oops::mpi::Scope _mpiScope(parentComm_);

  ASSERT_MSG(sourceParent.functionspace().mpi_comm() == parentComm_.name(),
    "CommStraightRedistribution::broadcastToSubMembers: Parent field does not "
    "have the same MPI communicator as the parent MPI communicator this redistribution has "
    "been set up with. Has the field been set up correctly?");

  ASSERT_MSG(targetSub.functionspace().mpi_comm() == subComm_.name(),
    "CommStraightRedistribution::broadcastToSubMembers: Subcomm field does not "
    "have the same MPI communicator as the sub MPI communicator this redistribution has "
    "been set up with. Has the field been set up correctly?");

  // Execute redistribution to the subcomm distributed field.
  const size_t numColors = colorInfo_.numColors();
  for (size_t color = 0; color < numColors; ++color) {
    transferFieldOnColor(sourceParent, color, targetSub);
  }

  oops::Log::trace() << "CommStraightRedistribution::broadcastToSubMembers done, from "
                     << parentComm_.name() << " to " << subComm_.name()
                     << std::endl;
}

// ------------------------------------------------------------------------------------------------

}  // namespace util

