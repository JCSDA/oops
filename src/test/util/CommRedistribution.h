/*
 * (C) Crown Copyright 2025-2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string_view>
#include <array>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"
#include "atlas/field/for_each.h"
#include "atlas/functionspace.h"
#include "atlas/grid.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/mpi/Comm.h"
#include "eckit/testing/Test.h"

#include "oops/../test/TestEnvironment.h"
#include "oops/mpi/ColorInfo.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/FieldSetOperations.h"
#include "oops/util/FunctionSpaceHelpers.h"
#include "oops/util/redistribution/CommRedistribution.h"

namespace test {

constexpr std::string_view parentCommName = "parent";
constexpr std::string_view subCommName = "subcommunicator";

// MPI Color group setups to run each test through.
// The WORLD communicator is split into a "parent" comm,
// where the size of this parent comm is the length of this array.
// The parent comm is then split further depening on the colors defined
// in these arrays.
static const std::vector<std::vector<size_t>> colorSetups {
  {0, 1},
  {0, 0, 1}, {0, 1, 1}, {0, 1, 2},
  {0, 0, 0, 1}, {0, 0, 1, 1}, {0, 1, 1, 1},
  {0, 1, 1, 2}, {0, 1, 2, 2},
  {0, 1, 2, 3},

  {0, 0, 0, 0, 0, 1}, {0, 0, 0, 0, 1, 1}, {0, 0, 0, 1, 1, 1},
  {0, 0, 1, 1, 1, 1}, {0, 1, 1, 1, 1, 1},

  {0, 0, 0, 0, 1, 2}, {0, 0, 0, 1, 1, 2}, {0, 0, 1, 1, 2, 2},
  {0, 0, 1, 2, 3, 3}, {0, 1, 2, 2, 3, 3}, {0, 1, 2, 3, 3, 3},
  {0, 1, 2, 3, 4, 4},
  {0, 1, 2, 3, 4, 5},
};

// ------------------------------------------------------------------------------------------------

// Confirm that the global indices are the same as the functionspace.
// (for test should be the case)
// Purpose being that if the data is in the correct place then this should match.
// A factor is given when a field is multiplied by a color, to reveal cross-talk
// issues across colors.
bool validateAgainstGidx(const atlas::Field& f,
                         const size_t globalHorizontalSize,
                         const double factor = 1.0) {
  const auto& fspace = f.functionspace();

  const auto gidx = atlas::array::make_view<atlas::uidx_t, 1>(fspace.global_index());
  const auto ghost = atlas::array::make_view<atlas::idx_t, 1>(fspace.ghost());
  const auto fview = atlas::array::make_view<double, 2>(f);
  ASSERT(gidx.shape(0) == fview.shape(0));

  bool passed = true;
  for (atlas::idx_t ij = 0; ij < gidx.shape(0); ++ij) {
    if (ghost(ij) == 0) {
      for (atlas::idx_t k = 0; k < fview.shape(1); ++k) {
        const double gidxDouble = static_cast<double>(gidx(ij) + (k * globalHorizontalSize))
            * factor;

        if (fview(ij, k) != gidxDouble) {
          passed = false;
          oops::Log::debug() << "MISMATCH (" << ij << ", " << k << ") => "
                             << fview(ij, k) << " != " << gidxDouble << std::endl;
        }
      }
    }
  }

  return passed;
}


/// This function:
/// 1) Creates a field containing the global index of each point,
///    including the vertical indexing.
///
/// 2) Redistributes the field to a subcommunicator with fewer MPI ranks.
///
/// 3) Compares the values in the redistributed field to the local global indices to ensure
///    the data matches as expected.
///
/// 4) Multiples the distributed field by a unique factor depending on the MPI color group
///    (for simplicity this is x *= (color+1)).
///
/// 4) Performs the opposite operation of redistributing to the "parent" communicator,
///    and ensures that each value from each color once again matches the global indices,
///    with the factor applied depending on the color.
///
void testRedistribution(const eckit::LocalConfiguration & config,
                        const eckit::mpi::Comm& comm = eckit::mpi::comm("world")) {
  oops::Log::trace() << "testRedistribution starting" << std::endl;

  // Setup functionspace on global communicator
  atlas::Grid grid{};
  atlas::grid::Partitioner partitioner{};
  atlas::Mesh mesh{};
  atlas::FunctionSpace fspace{};
  atlas::FieldSet fieldset{};
  util::setupFunctionSpace(comm, config, grid, partitioner,
                           mesh, fspace, fieldset);

  // Variables
  const std::string v = "var1";
  std::vector<std::string> varnames({v});
  eckit::LocalConfiguration varmeta, variablesconf;
  varmeta.set("levels", 2);
  variablesconf.set(varnames[0], varmeta);
  oops::Variables vars(variablesconf, varnames);

  // Create random fields
  atlas::FieldSet fset = util::createFieldSet(fspace, vars);

  // Get global horizontal size for indexing purposes.
  // Owned size...
  const size_t globalHorizontalSize = [&]() {
    const auto ghost = atlas::array::make_view<atlas::idx_t, 1>(fspace.ghost());

    size_t ownedSize = 0;
    for (atlas::idx_t ij = 0; ij < ghost.shape(0); ++ij) {
      ownedSize += 1 - ghost(ij);
    }

    comm.allReduceInPlace(ownedSize, eckit::mpi::sum());
    return ownedSize;
  }();

  // Set field to be equal to the global indices
  {
    // WARNING: PointCloud does not have a global_index field.
    auto parent = atlas::array::make_view<double, 2>(fset[v]);
    const auto gidx = atlas::array::make_view<atlas::uidx_t, 1>(fspace.global_index());

    for (atlas::idx_t ij = 0; ij < gidx.shape(0); ++ij) {
      for (atlas::idx_t k = 0; k < parent.shape(1); ++k) {
        parent(ij, k) = static_cast<double>(gidx(ij) + (k * globalHorizontalSize));
      }
    }
  }

  const double norm = util::normFieldSet(fset, vars.variables(), comm);

  const eckit::mpi::Comm & subComm = eckit::mpi::comm(subCommName);

  // Create function spaces on sub-communicator
  atlas::FunctionSpace subFspace;
  atlas::Grid subGrid{};
  atlas::grid::Partitioner subPartitioner{};
  atlas::Mesh subMesh{};
  atlas::FieldSet subFieldset{};
  util::setupFunctionSpace(subComm, config, subGrid, subPartitioner,
                           subMesh, subFspace, subFieldset);

  ASSERT(subFspace.size() >= fspace.size());

  const std::string redistName = config.getString("redistribution");
  // Setup redistribution object.
  const std::unique_ptr<util::CommRedistribution> redist(
      util::CommRedistributionFactory::create(redistName,
                                              subComm, comm, subFspace, fspace));

  atlas::FieldSet fsetsub = util::createFieldSet(subFspace, vars);
  util::zeroFieldSet(fsetsub);

  // -------- To subcomm --------
  redist->broadcastToSubMembers(fset[v], fsetsub[v]);

  // Check norm hasn't changed
  const double subdp = util::normFieldSet(fsetsub, vars.variables(), subComm);
  EXPECT(oops::is_close(subdp, norm, 1e-12));

  // Perform a scrict validation of the field values against the global indices.
  ASSERT_MSG(validateAgainstGidx(fsetsub[v], globalHorizontalSize),
             "Differences found after redistribution to sub communicator.");

  // Do something unique on each color to reveal any cross-talk/overwriting issues.
  // In this case multiply by the color.
  const oops::mpi::ColorInfo colInfo(subComm, comm);
  const double factor = static_cast<double>(colInfo.color() + 1);
  atlas::field::for_each_value(atlas::execution::par_unseq, fsetsub[v],
                               [&](double& x) {
                                 x *= factor;
                               });

  // -------- To parent --------
  const std::vector<atlas::Field> colorFields = redist->gatherToParent(fsetsub[v]);
  ASSERT(colorFields.size() == colInfo.numColors());

  size_t col = 0;
  for (const atlas::Field& f : colorFields) {
    ASSERT_MSG(validateAgainstGidx(f, globalHorizontalSize, static_cast<double>(col + 1)),
               "Differences found for color {" + std::to_string(col) + "}.");
    ++col;
  }

  oops::Log::trace() << "testRedistribution done." << std::endl;
}

// ------------------------------------------------------------------------------------------------

/// \brief This function calls the test with various MPI color setups.
void testRedistributionAtDiffSizes(const eckit::LocalConfiguration & config,
                                   const eckit::mpi::Comm& comm = eckit::mpi::comm("world")) {
  for (const std::vector<size_t>& colorSetup : colorSetups) {
    const size_t parentSize = colorSetup.size();
    ASSERT(parentSize <= comm.size());

    const size_t parentColor = comm.rank() < parentSize ? 0 : 1;
    comm.split(parentColor, std::string(parentCommName));

    // Any ranks in parent color 1 are not included in the test.
    // This lets the test run with different parent distribution sizes.
    if (parentColor == 0) {
      const eckit::mpi::Comm& parentComm = eckit::mpi::comm(parentCommName);

      // Split parent communicator.
      parentComm.split(colorSetup[parentComm.rank()],
                       std::string(subCommName));

      const oops::mpi::ColorInfo colInfo(eckit::mpi::comm(subCommName), parentComm);

      oops::Log::info() << "Testing with " << parentComm.size() << " parent ranks, "
                           "split: " << colInfo << std::endl;

      testRedistribution(config, parentComm);

      eckit::mpi::setCommDefault(parentCommName);
      eckit::mpi::deleteComm(subCommName);
    }

    eckit::mpi::setCommDefault(comm.name());
    eckit::mpi::deleteComm(parentCommName);
  }
}
// ------------------------------------------------------------------------------------------------
CASE("util/CommRedistribution/StructuredColumns/GatherScatter") {
  // Functionspace configuration
  eckit::LocalConfiguration config;
  config.set("function space", "StructuredColumns");
  config.set("grid.type", "regular_lonlat");
  config.set("grid.N", 20);
  config.set("halo", 1);
  config.set("partitioner", "checkerboard");
  config.set("redistribution", "gather-scatter");

  testRedistributionAtDiffSizes(config);
};

CASE("util/CommRedistribution/StructuredColumns/Straight") {
  // Functionspace configuration
  eckit::LocalConfiguration config;
  config.set("function space", "StructuredColumns");
  config.set("grid.type", "regular_lonlat");
  config.set("grid.N", 20);
  config.set("halo", 1);
  config.set("partitioner", "checkerboard");
  config.set("redistribution", "straight");

  testRedistributionAtDiffSizes(config);
};

CASE("util/CommRedistribution/NodeColumns/GatherScatter") {
  // Functionspace configuration
  eckit::LocalConfiguration config;
  config.set("function space", "NodeColumns");
  config.set("grid.name", "CS-LFR-15");
  config.set("partitioner", "cubedsphere");
  config.set("halo", 1);
  config.set("redistribution", "gather-scatter");

  testRedistributionAtDiffSizes(config);
};

CASE("util/CommRedistribution/NodeColumns/Straight") {
  // Functionspace configuration
  eckit::LocalConfiguration config;
  config.set("function space", "NodeColumns");
  config.set("grid.name", "CS-LFR-15");
  config.set("partitioner", "cubedsphere");
  config.set("halo", 1);
  config.set("redistribution", "straight");

  testRedistributionAtDiffSizes(config);
};


class CommRedistribution: public oops::Test {
 public:
  explicit CommRedistribution(const eckit::mpi::Comm& comm = oops::mpi::world()) :
    oops::Test(comm)
  {
    eckit::mpi::setCommDefault(comm.name());
  }
 private:
  std::string testid() const override {return "test::CommRedistribution";}

  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test
