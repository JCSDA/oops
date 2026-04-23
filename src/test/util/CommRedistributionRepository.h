/*
 * (C) Crown Copyright 2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string_view>
#include <memory>
#include <string>

#include "atlas/field.h"
#include "atlas/field/for_each.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/testing/Test.h"

#include "oops/../test/TestEnvironment.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/redistribution/CommRedistributionRepository.h"

namespace test {

constexpr std::string_view subCommName = "subcommunicator";

// ------------------------------------------------------------------------------------------------

// Throw away unneeded items when building functionspace.
atlas::FunctionSpace makeFSpace(const eckit::Configuration& config,
                                const eckit::mpi::Comm& comm) {
  atlas::Grid grid{};
  atlas::grid::Partitioner partitioner{};
  atlas::Mesh mesh{};
  atlas::FunctionSpace fspace{};
  atlas::FieldSet fieldset{};
  util::setupFunctionSpace(comm, config, grid, partitioner,
                           mesh, fspace, fieldset);
  return fspace;
}

// ------------------------------------------------------------------------------------------------

CASE("util/CommRedistributionRepository") {
  // Functionspace configuration
  eckit::LocalConfiguration config;
  config.set("function space", "StructuredColumns");
  config.set("grid.type", "regular_lonlat");
  config.set("grid.N", 10);
  config.set("halo", 1);
  config.set("partitioner", "checkerboard");

  const eckit::mpi::Comm& comm = eckit::mpi::comm("world");
  const eckit::mpi::Comm& subComm = eckit::mpi::comm(subCommName);

  // Doesn't matter much for the purposes of this test.
  const std::string redistName = "gather-scatter";

  // Setup functionspace on global communicator
  const atlas::FunctionSpace fspace = makeFSpace(config, comm);
  // Create functionspace on sub-communicator
  const atlas::FunctionSpace subFspace = makeFSpace(config, subComm);

  ASSERT(subFspace.size() >= fspace.size());

  // Request from the repository. This is the first time so it should be constructed.
  const util::CommRedistribution& redist =
      util::CommRedistributionRepository::get(redistName,
                                              subComm, comm, subFspace, fspace);

  // Do it again. This time, it should be fetched from the repository rather than creating a new
  // one. This can be checked by comparing pointer values...
  const util::CommRedistribution& redistFetched =
      util::CommRedistributionRepository::get(redistName,
                                              subComm, comm, subFspace, fspace);

  ASSERT(&redistFetched == &redist);


  // Now try with slightly different functionspace config...
  // Should expect the redistribution objects to be different.
  // Checks the key generation is unique enough.
  eckit::LocalConfiguration configEqualBands;
  configEqualBands.set("function space", "StructuredColumns");
  configEqualBands.set("grid.type", "regular_lonlat");
  configEqualBands.set("grid.N", 10);
  configEqualBands.set("halo", 1);
  configEqualBands.set("partitioner", "equal_bands");  // Different partitioner.

  const atlas::FunctionSpace fspaceEqualBands = makeFSpace(configEqualBands, comm);
  const atlas::FunctionSpace subFspaceEqualBands = makeFSpace(configEqualBands, subComm);

  const util::CommRedistribution& redistEqualBands =
      util::CommRedistributionRepository::get(redistName,
                                              subComm, comm,
                                              subFspaceEqualBands, fspaceEqualBands);

  ASSERT(&redistEqualBands != &redist);

  // Diff resolution
  eckit::LocalConfiguration configLowRes;
  configLowRes.set("function space", "StructuredColumns");
  configLowRes.set("grid.type", "regular_lonlat");
  configLowRes.set("grid.N", 5);
  configLowRes.set("halo", 1);
  configLowRes.set("partitioner", "checkerboard");

  const atlas::FunctionSpace fspaceLowRes = makeFSpace(configLowRes, comm);
  const atlas::FunctionSpace subFspaceLowRes = makeFSpace(configLowRes, subComm);

  const util::CommRedistribution& redistLowRes =
      util::CommRedistributionRepository::get(redistName,
                                              subComm, comm, subFspaceLowRes, fspaceLowRes);

  ASSERT(&redistLowRes != &redist);
}

// ------------------------------------------------------------------------------------------------

class CommRedistributionRepository: public oops::Test {
 public:
  explicit CommRedistributionRepository(const eckit::mpi::Comm& comm = oops::mpi::world()) :
    oops::Test(comm)
  {
    eckit::mpi::setCommDefault(comm.name());
    // Split comm in half.
    ASSERT(comm.size() % 2 == 0);
    comm.split(comm.rank() % 2, std::string(subCommName));
  }
 private:
  std::string testid() const override {return "test::CommRedistributionRepository";}

  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test
