/*
 * (C) Crown Copyright 2024-2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string_view>
#include <memory>
#include <random>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/testing/Test.h"

#include "oops/../test/TestEnvironment.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/FieldSetSubCommunicators.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/redistribution/CommRedistribution.h"

namespace test {

constexpr std::string_view subCommName = "subcommunicator";


size_t getColor(const eckit::mpi::Comm & comm) {
  return comm.rank() % 2;
}


void testSubCommunicators(const eckit::mpi::Comm & comm,
                          const eckit::LocalConfiguration & config) {
  oops::Log::trace() << "testSubCommunicators starting" << std::endl;
  // Assert that this communicator can be divided in two.
  ASSERT(comm.size() % 2 == 0);

  // Setup functionspace on global communicator
  atlas::Grid grid{};
  atlas::grid::Partitioner partitioner{};
  atlas::Mesh mesh{};
  atlas::FunctionSpace fspace{};
  atlas::FieldSet fieldset{};
  util::setupFunctionSpace(comm, config, grid, partitioner,
                           mesh, fspace, fieldset);

  // Variables
  eckit::LocalConfiguration variablesconf;
  std::vector<std::string> varnames;

  constexpr size_t numVars = 5;
  for (size_t ij = 0; ij < numVars; ++ij) {
    eckit::LocalConfiguration levelConf;
    levelConf.set("levels", (ij % 2)*2 + 2);  // alternate 2 and 4 levels.
    varnames.emplace_back("var" + std::to_string(ij));

    variablesconf.set(varnames.back(), levelConf);
  }
  oops::Variables vars(variablesconf, varnames);

  // Create random fields
  atlas::FieldSet fset = util::createRandomFieldSet(comm, fspace, vars);
  const double norm = util::normFieldSet(fset, vars.variables(), comm);

  comm.split(getColor(comm), std::string(subCommName));
  const eckit::mpi::Comm & subComm = eckit::mpi::comm(subCommName);

  // Create function spaces on sub-communicator
  atlas::FunctionSpace subFspace;
  atlas::Grid subGrid{};
  atlas::grid::Partitioner subPartitioner{};
  atlas::Mesh subMesh{};
  atlas::FieldSet subFieldset{};
  util::setupFunctionSpace(subComm, config, subGrid, subPartitioner,
                           subMesh, subFspace, subFieldset);

  // Create gather/scatter redistribution.
  const std::string redistName = config.getString("redistribution");
  // Copy fieldsets unto sub-communicators
  atlas::FieldSet fsetsub;
  util::redistributeToSubcommunicator(redistName, fset, fsetsub, subFspace);

  // Check norm hasn't changed
  const double subdp = util::normFieldSet(fsetsub, vars.variables(), subComm);
  EXPECT(oops::is_close(subdp, norm, 1e-12));

  // Do something unique on color 0.
  if (getColor(comm) == 0) {
    // Double
    util::multiplyFieldSet(fsetsub, 2.0);

    // Reorder the fields..!
    // There was an issue where one color changed the order of its fieldset
    // and it caused data to be mixed between the fields. This test will
    // pick this up.
    // Random order.
    std::random_device rd;
    std::mt19937 g(rd());

    std::vector<std::string> varsShuffled = fsetsub.field_names();
    std::shuffle(varsShuffled.begin(), varsShuffled.end(), g);

    atlas::FieldSet newFset;
    for (const std::string& fname : varsShuffled) {
      newFset.add(fsetsub[fname]);
    }
    // Replace fieldset
    fsetsub.clear();
    for (auto f : newFset) { fsetsub.add(f); }
  }

  // Gather and sum
  atlas::FieldSet fsetGathered;
  util::gatherAndSumFromSubcommunicator(redistName,
                                        fsetsub, fsetGathered,
                                        subFspace, fspace);

  // Compare to sum of initial fieldsets
  atlas::FieldSet doubleFset = util::copyFieldSet(fset);
  util::addFieldSets(doubleFset, fset);
  // NOTE(JC): Tripled.
  util::addFieldSets(doubleFset, fset);

  EXPECT(util::compareFieldSets(comm, fsetGathered, doubleFset));

  eckit::mpi::setCommDefault(comm.name());
  eckit::mpi::deleteComm(subComm.name());

  oops::Log::trace() << "testSubCommunicators done." << std::endl;
}

// ------------------------------------------------------------------------------------------------
CASE("util/FieldSetSubCommunicators/StructuredColumns/GatherScatter") {
  // Communicator
  const eckit::mpi::Comm & comm = oops::mpi::world();

  // Functionspace configuration
  eckit::LocalConfiguration config;
  config.set("function space", "StructuredColumns");
  config.set("grid.type", "regular_lonlat");
  config.set("grid.N", 20);
  config.set("halo", 1);
  config.set("redistribution", "gather-scatter");

  testSubCommunicators(comm, config);
};

CASE("util/FieldSetSubCommunicators/StructuredColumns/Straight") {
  // Communicator
  const eckit::mpi::Comm & comm = oops::mpi::world();

  // Functionspace configuration
  eckit::LocalConfiguration config;
  config.set("function space", "StructuredColumns");
  config.set("grid.type", "regular_lonlat");
  config.set("grid.N", 20);
  config.set("halo", 1);
  config.set("redistribution", "straight");

  testSubCommunicators(comm, config);
};

CASE("util/FieldSetSubCommunicators/NodeColumns/GatherScatter") {
  // Communicator
  const eckit::mpi::Comm & comm = oops::mpi::world();

  // Functionspace configuration
  eckit::LocalConfiguration config;
  config.set("function space", "NodeColumns");
  config.set("grid.name", "CS-LFR-15");
  config.set("partitioner", "cubedsphere");
  config.set("halo", 1);
  config.set("redistribution", "gather-scatter");

  testSubCommunicators(comm, config);
};

CASE("util/FieldSetSubCommunicators/NodeColumns/Straight") {
  // Communicator
  const eckit::mpi::Comm & comm = oops::mpi::world();

  // Functionspace configuration
  eckit::LocalConfiguration config;
  config.set("function space", "NodeColumns");
  config.set("grid.name", "CS-LFR-15");
  config.set("partitioner", "cubedsphere");
  config.set("halo", 1);
  config.set("redistribution", "straight");

  testSubCommunicators(comm, config);
};

class FieldSetSubCommunicators: public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "test::FieldSetSubCommunicators";}

  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test
