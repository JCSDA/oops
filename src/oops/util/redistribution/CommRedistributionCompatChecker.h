/*
 * (C) Crown Copyright 2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "eckit/mpi/Comm.h"


namespace util {

/// \brief Basic checks that can be shared for each CommRedistribution.
/// \details Includes checks such as comparing a lonlat hash between a given field
///          and ones the class has been set up with. Also included is a global
///          index check.
class CommRedistributionCompatChecker {
 public:
  CommRedistributionCompatChecker(const eckit::mpi::Comm& subComm,
                                  const eckit::mpi::Comm& parentComm,
                                  const atlas::FunctionSpace& subFSpace,
                                  const atlas::FunctionSpace& parentFSpace);

  void assertFSpaceMatchesSubSpace(const atlas::FunctionSpace& fspace) const;
  void assertFSpaceMatchesParentSpace(const atlas::FunctionSpace& fspace) const;

  void assertFieldMatchesSubSpace(const atlas::Field& f) const {
    return assertFSpaceMatchesSubSpace(f.functionspace());
  }
  void assertFieldMatchesParentSpace(const atlas::Field& f) const {
    return assertFSpaceMatchesParentSpace(f.functionspace());
  }

 private:
  const std::string subCommName_;
  const std::string parentCommName_;
  const atlas::FunctionSpace subFSpace_;
  const atlas::FunctionSpace parentFSpace_;

  const size_t subFSpaceOwnedSize_;
  const size_t parentFSpaceOwnedSize_;

  // Hashes for checking consistency
  const std::string subFSpaceLonLatHash_;
  const std::string parentFSpaceLonLatHash_;
};

}  // namespace util
