/*
 * (C) Copyright 2018-2019 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 */

#include "oops/util/IntSetParser.h"

#include <set>
#include <vector>

#include "eckit/testing/Test.h"

#include "oops/util/Logger.h"

namespace {

// -----------------------------------------------------------------------------

  CASE("test_get_channels_empty") {
    std::string chlist = "";
    std::set<int> expected{};

    std::set<int> channels = oops::parseIntSet(chlist);
    EXPECT(expected == channels);
  }

// -----------------------------------------------------------------------------

  CASE("test_get_one_channel") {
    std::string chlist = "5";
    std::set<int> expected{5};

    std::set<int> channels = oops::parseIntSet(chlist);
    EXPECT(expected == channels);
  }

// -----------------------------------------------------------------------------

  CASE("test_get_channel_list") {
    std::string chlist = "1,3,9";
    std::set<int> expected{1, 3, 9};

    std::set<int> channels = oops::parseIntSet(chlist);
    EXPECT(expected == channels);
  }

// -----------------------------------------------------------------------------

  CASE("test_get_identical_channels") {
    std::string chlist = "1,3,9,3";
    std::set<int> expected{1, 3, 9};

    std::set<int> channels = oops::parseIntSet(chlist);
    EXPECT(expected == channels);
  }

// -----------------------------------------------------------------------------

  CASE("test_get_channel_range") {
    std::string chlist = "3-7";
    std::set<int> expected{3, 4, 5, 6, 7};

    std::set<int> channels = oops::parseIntSet(chlist);
    EXPECT(expected == channels);
  }

// -----------------------------------------------------------------------------

  CASE("test_get_reverse_channel_range") {
    std::string chlist = "3-1";
    std::set<int> expected{};

    std::set<int> channels = oops::parseIntSet(chlist);
    EXPECT(expected == channels);
  }

// -----------------------------------------------------------------------------

  CASE("test_get_channels") {
    std::string chlist = "3-7,1,12-13";
    std::set<int> expected{1, 3, 4, 5, 6, 7, 12, 13};

    std::set<int> channels = oops::parseIntSet(chlist);
    EXPECT(expected == channels);
  }

}  // anonymous namespace

// -----------------------------------------------------------------------------

int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
