/*
 * (C) Crown Copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

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

namespace test {

void testSubCommunicators(const eckit::mpi::Comm & comm,
                          const eckit::LocalConfiguration & config) {
  oops::Log::trace() << "testSubCommunicators starting" << std::endl;

  // Setup functionspace on global communicator
  atlas::Grid grid{};
  atlas::grid::Partitioner partitioner{};
  atlas::Mesh mesh{};
  atlas::FunctionSpace fspace{};
  atlas::FieldSet fieldset{};
  util::setupFunctionSpace(comm, config, grid, partitioner,
                           mesh, fspace, fieldset);

  // Variables
  std::vector<std::string> varnames({"var1", "var2"});
  eckit::LocalConfiguration metaml2, metaml4, variablesconf;
  metaml2.set("levels", 2);
  metaml4.set("levels", 4);
  variablesconf.set(varnames[0], metaml2);
  variablesconf.set(varnames[1], metaml4);
  oops::Variables vars(variablesconf, varnames);

  // Create random fields
  atlas::FieldSet fset = util::createRandomFieldSet(comm, fspace, vars);
  const double norm = util::normFieldSet(fset, vars.variables(), comm);

  if (comm.size() == 2) {
    // Create split communicators
    const std::string commName = "subcommunicator";
    comm.split(comm.rank(), commName.c_str());
    const eckit::mpi::Comm & subComm = eckit::mpi::comm(commName.c_str());

    // Create function spaces on sub-communicator
    atlas::FunctionSpace subFspace;
    atlas::Grid grid{};
    atlas::grid::Partitioner partitioner{};
    atlas::Mesh mesh{};
    atlas::FieldSet fieldset{};
    util::setupFunctionSpace(subComm, config, grid, partitioner,
                             mesh, subFspace, fieldset);

    // Copy fieldsets unto sub-communicators
    atlas::FieldSet fsetsub;
    util::redistributeToSubcommunicator(fset, fsetsub,
                                        comm, subComm,
                                        fspace, subFspace);

    // Check norm hasn't changed
    const double subdp = util::normFieldSet(fsetsub, vars.variables(), subComm);
    EXPECT(oops::is_close(subdp, norm, 1e-12));

    // Gather and sum
    atlas::FieldSet fsetGathered;
    util::gatherAndSumFromSubcommunicator(fsetsub, fsetGathered,
                                          subComm, comm,
                                          subFspace, fspace);

    // Compare to sum of initial fieldsets
    atlas::FieldSet doubleFset = util::copyFieldSet(fset);
    util::addFieldSets(doubleFset, fset);
    EXPECT(util::compareFieldSets(comm, fsetGathered, doubleFset));

    eckit::mpi::setCommDefault(comm.name().c_str());
    eckit::mpi::deleteComm(subComm.name().c_str());
  }
  oops::Log::trace() << "testSubCommunicators done." << std::endl;
}

CASE("util/FieldSetSubCommunicators/StructuredColumns") {
  // Communicator
  const eckit::mpi::Comm & comm = oops::mpi::world();

  // Functionspace configuration
  eckit::LocalConfiguration config;
  config.set("function space", "StructuredColumns");
  config.set("grid.type", "regular_lonlat");
  config.set("grid.N", 10);
  config.set("halo", 1);

  testSubCommunicators(comm, config);
};

CASE("util/FieldSetSubCommunicators/NodeColumns") {
  // Communicator
  const eckit::mpi::Comm & comm = oops::mpi::world();

  // Functionspace configuration
  eckit::LocalConfiguration config;
  config.set("function space", "NodeColumns");
  config.set("grid.name", "CS-LFR-15");
  config.set("partitioner", "cubedsphere");
  config.set("halo", 1);

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
