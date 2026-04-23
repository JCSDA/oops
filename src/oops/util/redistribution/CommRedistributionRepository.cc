/*
 * (C) Crown Copyright 2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/Logger.h"
#include "oops/util/redistribution/CommRedistributionRepository.h"

namespace util {

// ------------------------------------------------------------------------------------------------
// Generate a unique id.
std::string CommRedistributionRepository::generateKey(const std::string& methodName,
                                                      const eckit::mpi::Comm& subComm,
                                                      const eckit::mpi::Comm& parentComm,
                                                      const atlas::FunctionSpace& subFSpace,
                                                      const atlas::FunctionSpace& parentFSpace) {
  const std::string sep = "_";
  return methodName + sep + subComm.name() + sep + parentComm.name() + sep +
         util::getGridUid(subFSpace) + sep +
         util::getLonLatHash(subFSpace) + util::getLonLatHash(parentFSpace);
}

// ------------------------------------------------------------------------------------------------

const CommRedistribution& CommRedistributionRepository::get(
    const std::string& methodName,
    const eckit::mpi::Comm& subComm,
    const eckit::mpi::Comm& parentComm,
    const atlas::FunctionSpace& subFSpace,
    const atlas::FunctionSpace& parentFSpace) {
  oops::Log::trace() << "CommRedistributionRepository::get starting." << std::endl;
  const std::string key = generateKey(methodName, subComm, parentComm, subFSpace, parentFSpace);

  auto& map = getInstance().instances_;
  const auto it = map.find(key);

  if (it != map.end()) {
    oops::Log::trace() << "CommRedistributionRepository::get finished." << std::endl;
    return *(it->second);
  } else {
    // Not in map - create a new one.
    map.emplace(key,
                CommRedistributionFactory::create(methodName,
                                                  subComm, parentComm,
                                                  subFSpace, parentFSpace));

    oops::Log::debug() << "CommRedistributionRepository::get: Created new entry for key '"
                       << key << "'." << std::endl;
    oops::Log::trace() << "CommRedistributionRepository::get finished." << std::endl;
    return *map.at(key);
  }
}

// ------------------------------------------------------------------------------------------------
}  // namespace util
