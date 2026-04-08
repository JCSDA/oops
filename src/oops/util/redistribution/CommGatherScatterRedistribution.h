/*
 * (C) Crown Copyright 2025-2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <vector>

#include "atlas/functionspace.h"

#include "oops/util/redistribution/CommRedistribution.h"
#include "oops/util/redistribution/CommRedistributionCompatChecker.h"

namespace util {

/// \brief Atlas Field 'gather-scatter' redistribution across communicators.
/// \details This redistribution works by gathering all data onto a single MPI rank,
///          and then scattering the data to the new functionspace, on another
///          MPI comm.
class CommGatherScatterRedistribution : public CommRedistribution {
 public:
  CommGatherScatterRedistribution(const eckit::mpi::Comm& sub,
                                  const eckit::mpi::Comm& parent,
                                  const atlas::FunctionSpace& subFSpace,
                                  const atlas::FunctionSpace& parentFSpace);

  void broadcastToSubMembers(const atlas::Field& sourceParent,
                             atlas::Field& targetSub) const override;

  std::vector<atlas::Field> gatherToParent(const atlas::Field& sourceSub) const override;

 private:
  static const std::string classname() { return "CommGatherScatterRedistribution"; }

  const eckit::mpi::Comm& subComm_;
  const eckit::mpi::Comm& parentComm_;

  const atlas::FunctionSpace subFSpace_;
  const atlas::FunctionSpace parentFSpace_;

  const CommRedistributionCompatChecker checker_;
};


}  // namespace util
