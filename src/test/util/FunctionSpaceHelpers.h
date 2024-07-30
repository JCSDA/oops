/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"
#include "atlas/grid/Partitioner.h"
#include "atlas/mesh.h"
#include "atlas/meshgenerator.h"

#include "oops/../test/TestEnvironment.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/missingValues.h"

namespace test {

// Test of StructuredMeshToStructuredColumnsIndexMap class
// Loop over every mesh point (including halo), and check it is mapped either to,
// - the same point (=> same lon,lat) in the function space, OR
// - to a missing value, denoting a mesh point that is not in the function space (on this MPI task)
void testIndexMapper(const atlas::Mesh & mesh, const atlas::FunctionSpace & fs) {
  const atlas::functionspace::StructuredColumns structuredcolumns(fs);
  if (structuredcolumns) {
    // check every mesh point is mapped to the same coordinate, or is not mapped
    auto fs_lonlat = atlas::array::make_view<double, 2>(structuredcolumns.lonlat());
    auto mesh_lonlat = atlas::array::make_view<double, 2>(mesh.nodes().lonlat());

    util::StructuredMeshToStructuredColumnsIndexMap map{};
    map.initialize(mesh, structuredcolumns);

    for (atlas::idx_t imesh = 0; imesh < mesh_lonlat.shape(0); ++imesh) {
      const atlas::idx_t ifs = map(imesh);
      if (ifs != util::missingValue<atlas::idx_t>()) {
        EXPECT(mesh_lonlat(imesh, 0) == fs_lonlat(ifs, 0));
        EXPECT(mesh_lonlat(imesh, 1) == fs_lonlat(ifs, 1));
      }
    }
  }
}

void testStructuredColumnsLonLat(const eckit::mpi::Comm & comm) {
  eckit::LocalConfiguration config;
  config.set("function space", "StructuredColumns");
  config.set("grid.type", "regular_lonlat");
  config.set("grid.N", 10);
  config.set("halo", 1);

  atlas::Grid grid{};
  atlas::grid::Partitioner partitioner{};
  atlas::Mesh mesh{};
  atlas::FunctionSpace functionspace{};
  atlas::FieldSet fieldset{};

  util::setupFunctionSpace(comm, config, grid, partitioner, mesh, functionspace, fieldset);

  EXPECT(grid.type() == "structured");
  EXPECT(grid.name() == "L40x21");
  EXPECT(grid.uid() == "53ed3a5551fa8cdb079fe115519e0a92");
  EXPECT(partitioner.type() == "equal_regions");
  EXPECT(functionspace.type() == "StructuredColumns");
  EXPECT(fieldset.has("owned"));

  testIndexMapper(mesh, functionspace);
}

CASE("util/FunctionSpaceHelpers/StructuredColumnsLonLat") {
  const eckit::mpi::Comm & comm = oops::mpi::world();
  testStructuredColumnsLonLat(comm);

  if (comm.size() > 1) {
    const std::string commName = "subcommunicator";
    comm.split(comm.rank(), commName.c_str());
    testStructuredColumnsLonLat(eckit::mpi::comm(commName.c_str()));
    eckit::mpi::setCommDefault(comm.name().c_str());
    eckit::mpi::deleteComm(commName.c_str());
  }
}

void testStructuredColumnsGaussian(const eckit::mpi::Comm & comm) {
  eckit::LocalConfiguration config;
  config.set("function space", "StructuredColumns");
  config.set("grid.type", "regular_gaussian");
  config.set("grid.N", 10);
  config.set("halo", 1);

  atlas::Grid grid{};
  atlas::grid::Partitioner partitioner{};
  atlas::Mesh mesh{};
  atlas::FunctionSpace functionspace{};
  atlas::FieldSet fieldset{};

  util::setupFunctionSpace(comm, config, grid, partitioner, mesh, functionspace, fieldset);

  EXPECT(grid.type() == "structured");
  EXPECT(grid.name() == "F10");
  EXPECT(grid.uid() == "2734f1f878e2e047d290b3a578fc2927");
  EXPECT(partitioner.type() == "equal_regions");
  EXPECT(functionspace.type() == "StructuredColumns");
  EXPECT(fieldset.has("owned"));

  testIndexMapper(mesh, functionspace);
}

CASE("util/FunctionSpaceHelpers/StructuredColumnsGaussian") {
  const eckit::mpi::Comm & comm = oops::mpi::world();
  testStructuredColumnsGaussian(comm);

  if (comm.size() > 1) {
    const std::string commName = "subcommunicator";
    comm.split(comm.rank(), commName.c_str());
    testStructuredColumnsGaussian(eckit::mpi::comm(commName.c_str()));
    eckit::mpi::setCommDefault(comm.name().c_str());
    eckit::mpi::deleteComm(commName.c_str());
  }
}

void testStructuredColumnsGaussianCustomDistribution(const eckit::mpi::Comm & comm) {
  eckit::LocalConfiguration config;
  config.set("function space", "StructuredColumns");
  config.set("grid.type", "regular_gaussian");
  config.set("grid.N", 10);
  config.set("halo", 0);
  config.set("no point on last task", true);

  atlas::Grid grid{};
  atlas::grid::Partitioner partitioner{};
  atlas::Mesh mesh{};
  atlas::FunctionSpace functionspace{};
  atlas::FieldSet fieldset{};

  util::setupFunctionSpace(comm, config, grid, partitioner, mesh, functionspace, fieldset);

  EXPECT(grid.type() == "structured");
  EXPECT(grid.name() == "F10");
  EXPECT(grid.uid() == "2734f1f878e2e047d290b3a578fc2927");
  EXPECT(partitioner.type() == "equal_regions");
  EXPECT(partitioner.nb_partitions() == std::max(2ul, comm.size()) - 1);
  EXPECT(functionspace.type() == "StructuredColumns");
  EXPECT(fieldset.has("owned"));

  testIndexMapper(mesh, functionspace);
}

CASE("util/FunctionSpaceHelpers/StructuredColumnsGaussianCustomDistribution") {
  const eckit::mpi::Comm & comm = oops::mpi::world();
  testStructuredColumnsGaussianCustomDistribution(comm);

  if (comm.size() > 1) {
    const std::string commName = "subcommunicator";
    comm.split(comm.rank(), commName.c_str());
    testStructuredColumnsGaussianCustomDistribution(eckit::mpi::comm(commName.c_str()));
    eckit::mpi::setCommDefault(comm.name().c_str());
    eckit::mpi::deleteComm(commName.c_str());
  }
}

void testStructuredColumnsRegional(const eckit::mpi::Comm & comm) {
  eckit::LocalConfiguration config;
  config.set("function space", "StructuredColumns");
  config.set("grid.type", "regional");
  config.set("grid.nx", 71);
  config.set("grid.ny", 53);
  config.set("grid.dx", 2.5e3);
  config.set("grid.dy", 2.5e3);
  config.set("grid.lonlat(centre)", std::vector<double>{{9.9, 56.3}});
  config.set("grid.projection.type", "lambert_conformal_conic");
  config.set("grid.projection.longitude0", 0.0);
  config.set("grid.projection.latitude0", 56.3);
  config.set("grid.y_numbering", 1);
  config.set("partitioner", "checkerboard");

  atlas::Grid grid{};
  atlas::grid::Partitioner partitioner{};
  atlas::Mesh mesh{};
  atlas::FunctionSpace functionspace{};
  atlas::FieldSet fieldset{};

  util::setupFunctionSpace(comm, config, grid, partitioner, mesh, functionspace, fieldset);

  EXPECT(grid.type() == "structured");  // NB: not "regional" !
  EXPECT(grid.name() == "structured");
  EXPECT(grid.uid() == "60797064f97ff3149d4814f6105028d5");
  EXPECT(partitioner.type() == "checkerboard");
  EXPECT(functionspace.type() == "StructuredColumns");
  EXPECT(fieldset.has("owned"));

  testIndexMapper(mesh, functionspace);
}

CASE("util/FunctionSpaceHelpers/StructuredColumnsRegional") {
  const eckit::mpi::Comm & comm = oops::mpi::world();
  testStructuredColumnsRegional(comm);

  if (comm.size() > 1) {
    const std::string commName = "subcommunicator";
    comm.split(comm.rank(), commName.c_str());
    testStructuredColumnsRegional(eckit::mpi::comm(commName.c_str()));
    eckit::mpi::setCommDefault(comm.name().c_str());
    eckit::mpi::deleteComm(commName.c_str());
  }
}

void testNodeColumnsCubedSphere(const eckit::mpi::Comm & comm) {
  eckit::LocalConfiguration config;
  config.set("function space", "NodeColumns");
  config.set("grid.name", "CS-LFR-12");
  config.set("partitioner", "cubedsphere");
  config.set("halo", 1);

  atlas::Grid grid{};
  atlas::grid::Partitioner partitioner{};
  atlas::Mesh mesh{};
  atlas::FunctionSpace functionspace{};
  atlas::FieldSet fieldset{};

  util::setupFunctionSpace(comm, config, grid, partitioner, mesh, functionspace, fieldset);

  EXPECT(grid.type() == "cubedsphere");
  EXPECT(grid.name() == "CS-LFR-12");
  EXPECT(grid.uid() == "8aa0b472107ce06c53b5c760886b9fb1");
  EXPECT(partitioner.type() == "cubedsphere");
  EXPECT(functionspace.type() == "NodeColumns");
  EXPECT(fieldset.has("owned"));
}

CASE("util/FunctionSpaceHelpers/NodeColumnsCubedSphere") {
  const eckit::mpi::Comm & comm = oops::mpi::world();
  testNodeColumnsCubedSphere(comm);

  if (comm.size() > 1) {
    const std::string commName = "subcommunicator";
    comm.split(comm.rank(), commName.c_str());
    testNodeColumnsCubedSphere(eckit::mpi::comm(commName.c_str()));
    eckit::mpi::setCommDefault(comm.name().c_str());
    eckit::mpi::deleteComm(commName.c_str());
  }
}

void testNodeColumnsUnstructured(const eckit::mpi::Comm & comm) {
  eckit::LocalConfiguration config;
  config.set("function space", "NodeColumns");
  config.set("grid.type", "unstructured");
  config.set("grid.xy", std::vector<double>{{0, 0, 90, 0, 180, 0, 270, 0, 0, 90, 0, -90}});
  config.set("partitioner", "equal_regions");
  config.set("no point on last task", true);

  atlas::Grid grid{};
  atlas::grid::Partitioner partitioner{};
  atlas::Mesh mesh{};
  atlas::FunctionSpace functionspace{};
  atlas::FieldSet fieldset{};

  util::setupFunctionSpace(comm, config, grid, partitioner, mesh, functionspace, fieldset);

  EXPECT(grid.type() == "unstructured");
  EXPECT(grid.name() == "unstructured.7289032");
  EXPECT(grid.uid() == "7289032475be87656a165901243eda16");
  EXPECT(partitioner.type() == "equal_regions");
  EXPECT(functionspace.type() == "NodeColumns");
  EXPECT(fieldset.has("owned"));
}

CASE("util/FunctionSpaceHelpers/NodeColumnsUnstructured") {
  const eckit::mpi::Comm & comm = oops::mpi::world();
  testNodeColumnsUnstructured(comm);

  if (comm.size() > 1) {
    const std::string commName = "subcommunicator";
    comm.split(comm.rank(), commName.c_str());
    testNodeColumnsUnstructured(eckit::mpi::comm(commName.c_str()));
    eckit::mpi::setCommDefault(comm.name().c_str());
    eckit::mpi::deleteComm(commName.c_str());
  }
}

class FunctionSpaceHelpers : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "test::FunctionSpaceHelpers";}

  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test
