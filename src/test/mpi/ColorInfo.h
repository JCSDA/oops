/*
 * (C) Crown Copyright 2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string_view>
#include <array>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
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

  // ===== Unconventional setups. Not contiguous. =====
  {0, 1, 0},
  {0, 0, 1, 0},
  {0, 1, 0, 1},
  {0, 1, 2, 0},
  {0, 1, 2, 1},
  {0, 1, 1, 2, 0},
  {0, 1, 2, 3, 0},
  {0, 1, 0, 1, 0, 1},
  {0, 1, 2, 1, 0, 1},
  {0, 1, 2, 1, 3, 0},

  // ===== Setups that ignore assumptions. NOT ADVISED. =====
  // The problem with these setups is that they introduce new colors out of order, or not
  // starting from 0.
  //
  // This information is lost when trying to infer the colors because as colors are discovered
  // by `ColorInfo` they are given a new index. This is perfectly fine in the cases when
  // colors are introduced in order (as is normal).
  //
  // This has not been made an error case in this test because this is expected output from
  // this object. There is no way for ColorInfo to infer the original keys used for splitting
  // from the communicators alone, so there would be no way to throw an error in this case.

  //               Output expected from ColorInfo
  {1, 0},          // {0, 1}
  {2, 1, 0},       // {0, 1, 2}
  {1, 0, 1},       // {0, 1, 0}
  {12, 63, 2},     // {0, 1, 2}
};

// ------------------------------------------------------------------------------------------------

std::vector<size_t> calculateExpectedRoots(const std::vector<size_t>& colorSetup) {
  std::vector<size_t> roots;
  std::set<size_t> colorsSeen;

  size_t parentRank = 0;
  for (const size_t col : colorSetup) {
    if (colorsSeen.find(col) == colorsSeen.end()) {
      roots.push_back(parentRank);
      colorsSeen.insert(col);
    }
    ++parentRank;
  }

  return roots;
}

// ------------------------------------------------------------------------------------------------

// Generates what ColorInfo will see the colors as when gathering information via
// the MPI communicators.
std::vector<size_t> calculateExpectedColorOutput(const std::vector<size_t>& colorSetup) {
  std::unordered_map<size_t, size_t> mapping;  // from color setup -> expected color.

  // Translated colors.
  std::vector<size_t> colors;
  colors.reserve(colorSetup.size());

  size_t colorId = 0;
  for (const size_t col : colorSetup) {
    // Unseen color?
    if (mapping.find(col) == mapping.end()) {
      mapping.emplace(col, colorId++);
    }

    colors.emplace_back(mapping[col]);
  }

  ASSERT(colors.size() == colorSetup.size());
  return colors;
}

// ------------------------------------------------------------------------------------------------

CASE("mpi/ColorInfo") {
  const eckit::mpi::Comm& comm = eckit::mpi::comm("world");

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

      oops::Log::info() << "----------------------------------------------------------------------"
                        << std::endl
                        << "- Running setup => " << colorSetup << std::endl;

      const std::vector<size_t> expectedRoots = calculateExpectedRoots(colorSetup);
      const std::vector<size_t> expectedColors = calculateExpectedColorOutput(colorSetup);
      const oops::mpi::ColorInfo colInfo(subCommName, parentCommName);

      oops::Log::info() << "- " << colInfo << std::endl;

      oops::Log::info() << "- Roots => " << colInfo.roots() << " expected "
                        << expectedRoots << std::endl
                        << "- Colour => " << colInfo.color() << " expected "
                        << colorSetup[parentComm.rank()] << std::endl;
      EXPECT(expectedRoots == colInfo.roots());
      EXPECT(expectedRoots.size() == colInfo.numColors());
      EXPECT(expectedColors[parentComm.rank()] == colInfo.color());

      eckit::mpi::setCommDefault(parentCommName);
      eckit::mpi::deleteComm(subCommName);
    }

    eckit::mpi::setCommDefault(comm.name());
    eckit::mpi::deleteComm(parentCommName);
  }
};

// ------------------------------------------------------------------------------------------------

class ColorInfo: public oops::Test {
 public:
  explicit ColorInfo(const eckit::mpi::Comm& comm = oops::mpi::world()) :
    oops::Test(comm)
  {
    eckit::mpi::setCommDefault(comm.name());
  }
 private:
  std::string testid() const override {return "test::ColorInfo";}

  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test
