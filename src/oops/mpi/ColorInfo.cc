/*
 * (C) Crown Copyright 2025-2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <algorithm>

#include "eckit/exception/Exceptions.h"

#include "oops/mpi/ColorInfo.h"
#include "oops/util/Logger.h"
#include "oops/util/Timer.h"

namespace oops {
namespace mpi {

ColorInfo::ColorInfo(const eckit::mpi::Comm& subComm,
                     const eckit::mpi::Comm& parentComm) :
  parentRoots_(parentComm.size()),
  parentName_(parentComm.name()),
  subName_(subComm.name())
{
  const util::Timer _timer(classname(), classname());
  Log::trace() << "ColorInfo::ColorInfo starting." << std::endl;

  constexpr size_t subRoot = 0;

  ASSERT(subComm.size() < parentComm.size());

  // Rank of local root on parent communicator
  size_t rootOnParent = 0;
  if (subComm.rank() == subRoot) {
    rootOnParent = parentComm.rank();
  }

  // Share parent root value within sub-communicator
  subComm.broadcast(rootOnParent, subRoot);

  // Gather local root ranks from all sub-communicators
  parentComm.allGather(rootOnParent, parentRoots_.begin(), parentRoots_.end());

  // Only keep one instance of each local root
  std::sort(parentRoots_.begin(), parentRoots_.end());
  auto last = std::unique(parentRoots_.begin(), parentRoots_.end());
  parentRoots_.erase(last, parentRoots_.end());

  // Get index of current sub-communicator in vector of `roots`
  color_ = std::distance(parentRoots_.begin(),
                         std::find(parentRoots_.begin(), parentRoots_.end(), rootOnParent));

  Log::trace() << "ColorInfo::ColorInfo finished." << std::endl;
}

// ------------------------------------------------------------------------------------------------

ColorInfo::ColorInfo(const std::string_view subCommName,
                     const std::string_view parentCommName) :
  ColorInfo(eckit::mpi::comm(subCommName),
             eckit::mpi::comm(parentCommName))
{}

// ------------------------------------------------------------------------------------------------

void ColorInfo::print(std::ostream & os) const {
  os << "ColorInfo(" << subName_ << ", " << parentName_ << ") { parent rank: "
     << eckit::mpi::comm(parentName_).rank()
     << ", colors: " << numColors()
     << ", this: " << color()
     << " }";
}

// ------------------------------------------------------------------------------------------------

}  // namespace mpi
}  // namespace oops

