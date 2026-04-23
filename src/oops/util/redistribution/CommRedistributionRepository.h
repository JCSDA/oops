/*
 * (C) Crown Copyright 2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>
#include <unordered_map>

#include "atlas/functionspace.h"
#include "eckit/mpi/Comm.h"
#include "oops/util/redistribution/CommRedistribution.h"

namespace util {

/// \brief Multiton repository holding `CommRedistribution` objects.
///        A user requests a given redistribution object. If it exists,
///        it is returned. If it does not already exist, a new one is created
///        and added to the repository for later use.
class CommRedistributionRepository {
 public:
  static const CommRedistribution& get(const std::string& methodName,
                                       const eckit::mpi::Comm& subComm,
                                       const eckit::mpi::Comm& parentComm,
                                       const atlas::FunctionSpace& subFSpace,
                                       const atlas::FunctionSpace& parentFSpace);

 private:
  CommRedistributionRepository() {}
  CommRedistributionRepository(const CommRedistributionRepository&)            = delete;
  CommRedistributionRepository& operator=(const CommRedistributionRepository&) = delete;

  static std::string generateKey(const std::string& methodName,
                                 const eckit::mpi::Comm& subComm,
                                 const eckit::mpi::Comm& parentComm,
                                 const atlas::FunctionSpace& subFSpace,
                                 const atlas::FunctionSpace& parentFSpace);

  static CommRedistributionRepository& getInstance() {
    static CommRedistributionRepository theInstance;
    return theInstance;
  }

  std::unordered_map<std::string, std::unique_ptr<CommRedistribution>> instances_;
};

}  // namespace util
