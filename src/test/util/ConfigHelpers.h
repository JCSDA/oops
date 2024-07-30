/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/runs/Test.h"
#include "oops/util/ConfigHelpers.h"
#include "oops/util/Expect.h"

namespace test {

CASE("util/ConfigHelpers/setMember") {
  // Check setMember without member template
  eckit::LocalConfiguration conf{};
  util::setMember(conf, 7);
  EXPECT(conf.getInt("member") == 7);

  // Check setMember with member template
  conf.set("member pattern", "%mem%");
  conf.set("filename", "path/to/foo%mem%.txt");
  util::setMember(conf, 3);
  EXPECT(conf.getString("filename") == "path/to/foo3.txt");

  EXPECT(conf.getInt("member") == 7);  // Note this is unchanged
}

class ConfigHelpers : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "test::ConfigHelpers";}

  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test
