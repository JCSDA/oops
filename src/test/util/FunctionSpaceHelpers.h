/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

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

namespace test {

CASE("util/FunctionSpaceHelpers/StructuredColumnsLonLat") {
  const eckit::mpi::Comm & comm = oops::mpi::world();

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
}

CASE("util/FunctionSpaceHelpers/StructuredColumnsGaussian") {
  const eckit::mpi::Comm & comm = oops::mpi::world();

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
}

CASE("util/FunctionSpaceHelpers/StructuredColumnsRegional") {
  const eckit::mpi::Comm & comm = oops::mpi::world();

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

  std::cout << grid.type() << std::endl;
  std::cout << grid.name() << std::endl;
  std::cout << grid.uid() << std::endl;
  std::cout << partitioner.type() << std::endl;
  std::cout << functionspace.type() << std::endl;
  EXPECT(grid.type() == "structured");  // NB: not "regional" !
  EXPECT(grid.name() == "structured");
  EXPECT(grid.uid() == "60797064f97ff3149d4814f6105028d5");
  EXPECT(partitioner.type() == "checkerboard");
  EXPECT(functionspace.type() == "StructuredColumns");
  EXPECT(fieldset.has("owned"));
}

CASE("util/FunctionSpaceHelpers/NodeColumnsCubedSphere") {
  const eckit::mpi::Comm & comm = oops::mpi::world();

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

  std::cout << grid.type() << std::endl;
  std::cout << grid.name() << std::endl;
  std::cout << grid.uid() << std::endl;
  std::cout << partitioner.type() << std::endl;
  std::cout << functionspace.type() << std::endl;
  EXPECT(grid.type() == "cubedsphere");
  EXPECT(grid.name() == "CS-LFR-12");
  EXPECT(grid.uid() == "8aa0b472107ce06c53b5c760886b9fb1");
  EXPECT(partitioner.type() == "cubedsphere");
  EXPECT(functionspace.type() == "NodeColumns");
  EXPECT(fieldset.has("owned"));
}

CASE("util/FunctionSpaceHelpers/NodeColumnsUnstructured") {
  const eckit::mpi::Comm & comm = oops::mpi::world();

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

class FunctionSpaceHelpers : public oops::Test {
 private:
  std::string testid() const override {return "test::FunctionSpaceHelpers";}

  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test
