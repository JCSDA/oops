/*
 * (C) Crown Copyright 2025-2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "eckit/exception/Exceptions.h"

#include "oops/mpi/ColorInfo.h"
#include "oops/mpi/Scope.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/redistribution/CommGatherScatterRedistribution.h"
#include "oops/util/Timer.h"

namespace util {

// Register in factory.
static CommRedistributionFactory::Maker<CommGatherScatterRedistribution>
    makerGatherScatterRedist_("gather-scatter");

// -----------------------------------------------------------------------------

CommGatherScatterRedistribution::CommGatherScatterRedistribution(
    const eckit::mpi::Comm& sub,
    const eckit::mpi::Comm& parent,
    const atlas::FunctionSpace& subFSpace,
    const atlas::FunctionSpace& parentFSpace) :
  subComm_(sub), parentComm_(parent),
  subFSpace_(subFSpace), parentFSpace_(parentFSpace),
  checker_(sub, parent, subFSpace, parentFSpace)
{
  Timer _timer(classname(), "CommGatherScatterRedistribution");
}


// -----------------------------------------------------------------------------

void CommGatherScatterRedistribution::broadcastToSubMembers(const atlas::Field& sourceParent,
                                                atlas::Field& targetSub) const {
  Timer _timer(classname(), "broadcastToSubMembers");
  const oops::mpi::Scope _parentScope(parentComm_);

  oops::Log::trace() << "CommGatherScatterRedistribution::broadcastToSubMembers starting, from "
                     << parentComm_.name() << " to " << subComm_.name() << std::endl;

  checker_.assertFieldMatchesParentSpace(sourceParent);
  checker_.assertFieldMatchesSubSpace(targetSub);
  ASSERT(sourceParent.levels() == targetSub.levels());

  // 1. Gather fields onto global fields on each sub communicator, on parent comm.
  // TODO(Mayeul/JC): add different data types (int, long, float)
  ASSERT(sourceParent.datatype() == atlas::array::DataType::kind<double>());

  const auto fieldConfigBase = atlas::option::name(sourceParent.name()) |
                               atlas::option::levels(sourceParent.levels()) |
                               atlas::option::datatype(sourceParent.datatype());

  const atlas::FunctionSpace& fspaceIn = sourceParent.functionspace();
  const atlas::FunctionSpace& fspaceOut = targetSub.functionspace();

  const oops::mpi::ColorInfo colorInfo(subComm_, parentComm_);

  // Gather on each sub-communicator
  std::vector<atlas::Field> rootFields;
  for (const size_t otherRoot : colorInfo.roots()) {
    rootFields.push_back(fspaceIn.createField(fieldConfigBase |
                                              atlas::option::global(otherRoot)));

    fspaceIn.gather(sourceParent, rootFields.back());
  }

  // 2. Copy global fields defined on global communicator to global field
  //    defined on sub-communicator
  {
    const oops::mpi::Scope _subScope(subComm_);
    // Changing the default communicator means we should change the owner of the global field
    // from root to sub-root. Unfortunately, atlas field offers no way to change this.
    // As an alternative, we copy the field into another global field on the same PE
    // with the correct root numbering.
    atlas::Field rootField = fspaceOut.createField(fieldConfigBase |
                                                   atlas::option::global(0));
    atlas::array::make_view<double, 2>(rootField).assign(
                atlas::array::make_view<double, 2>(rootFields[colorInfo.color()]));

    // 3. Scatter from local root PEs to each subcommunicator
    fspaceOut.scatter(rootField, targetSub);

    targetSub.set_dirty();
  }

  oops::Log::trace() << "CommGatherScatterRedistribution::broadcastToSubMembers done, from "
                     << parentComm_.name() << " to " << subComm_.name()
                     << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<atlas::Field> CommGatherScatterRedistribution::gatherToParent(
    const atlas::Field& sourceSub) const {
  Timer _timer(classname(), "gatherToParent");
  const oops::mpi::Scope _parentScope(parentComm_);

  oops::Log::trace() << "CommGatherScatterRedistribution::gatherToParent starting, from "
                     << subComm_.name() << " to " << parentComm_.name()
                     << std::endl;

  checker_.assertFieldMatchesSubSpace(sourceSub);

  // TODO(Mayeul/JC): add different data types (int, long, float)
  ASSERT(sourceSub.datatype() == atlas::array::DataType::kind<double>());

  const auto fieldConfigBase = atlas::option::name(sourceSub.name()) |
                               atlas::option::levels(sourceSub.levels()) |
                               atlas::option::datatype(sourceSub.datatype());

  // On this current split comm, gather the field on the root of the split communicator.
  const atlas::Field subCommGlobalField = [&]() {
    const oops::mpi::Scope _subScope(subComm_);

    const atlas::FunctionSpace& fspaceIn = sourceSub.functionspace();
    atlas::Field subGField = fspaceIn.createField(fieldConfigBase |
                                                  atlas::option::global(0));

    fspaceIn.gather(sourceSub, subGField);
    return subGField;
  }();

  // Copy sub-communicator global field to global-communicator global field,
  // for each sub split communicator.
  // (back on parent communicator)
  const oops::mpi::ColorInfo colorInfo(subComm_, parentComm_);

  std::vector<atlas::Field> globallyDistributedFields{};
  for (const size_t otherRoot : colorInfo.roots()) {
    // Global field on the root of the split comm
    atlas::Field parentCommGlobalField = parentFSpace_.createField(
        fieldConfigBase |
        atlas::option::global(otherRoot));

    // Copy the global field on sub comm to the global field on parent comm.
    const auto subGFieldView = atlas::array::make_view<const double, 2>(subCommGlobalField);
    atlas::array::make_view<double, 2>(parentCommGlobalField).assign(subGFieldView);

    // Distribute the copied global field on the parent comm.
    atlas::Field parentField = parentFSpace_.createField(fieldConfigBase);

    parentFSpace_.scatter(parentCommGlobalField, parentField);

    globallyDistributedFields.push_back(parentField);
  }

  oops::Log::trace() << "CommGatherScatterRedistribution::gatherToParent done, from "
                     << subComm_.name() << " to " << parentComm_.name()
                     << std::endl;
  return globallyDistributedFields;
}

// -----------------------------------------------------------------------------
}  // namespace util

