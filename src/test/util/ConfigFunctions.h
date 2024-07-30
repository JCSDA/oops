/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_CONFIGFUNCTIONS_H_
#define TEST_UTIL_CONFIGFUNCTIONS_H_

#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/config/YAMLConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/runs/Test.h"
#include "oops/util/ConfigFunctions.h"
#include "oops/util/Expect.h"

namespace test {

  eckit::LocalConfiguration createConfig() {
    // Main configuration
    eckit::LocalConfiguration conf;

    // Subconfiguration
    std::string subStr("{\"int\": 1, \"real\": 1.0,\"string\": \"dummy_%pattern%\"}");
    const eckit::YAMLConfiguration yamlSub(subStr);
    const eckit::LocalConfiguration sub(yamlSub);
    conf.set("sub", sub);

    // Final pair
    eckit::LocalConfiguration fp(sub, "int");
    conf.set("fp", fp);

    // Vector configuration
    std::string vectorStr("[{\"int\": 1, \"real\": 1.0,\"string\": \"dummy_%pattern%\"}]");
    const eckit::YAMLConfiguration yamlVector(vectorStr);
    const eckit::LocalConfiguration vector(yamlVector);
    conf.set("vector", vector);

    // Reference subconfiguration with "%pattern%" replaced by "001"
    std::string refStr("{\"int\": 1, \"real\": 1.0,\"string\": \"dummy_001\"}");
    const eckit::YAMLConfiguration yamlRef(refStr);
    const eckit::LocalConfiguration ref(yamlRef);
    conf.set("ref", ref);

    // Other subconfiguration
    std::string otherStr("{\"other int\": 2}");
    const eckit::YAMLConfiguration yamlOther(otherStr);
    const eckit::LocalConfiguration other(yamlOther);
    conf.set("other", other);

    // Merged configuration
    std::string mergedStr("{\"int\": 1, \"real\": 1.0,\"string\": \"dummy_001\",\"other int\": 2}");
    const eckit::YAMLConfiguration yamlMerged(mergedStr);
    const eckit::LocalConfiguration merged(yamlMerged);
    conf.set("merged", merged);

    return conf;
  }

CASE("util/ConfigFunctions/isVector") {
  eckit::LocalConfiguration conf = createConfig();

  const eckit::LocalConfiguration sub(conf, "sub");
  EXPECT_NOT(util::isVector(sub));

  const eckit::LocalConfiguration fp(conf, "fp");
  EXPECT_NOT(util::isVector(fp));

  const eckit::LocalConfiguration vector(conf, "vector");
  EXPECT(util::isVector(vector));
}

CASE("util/ConfigFunctions/isSubConfig") {
  eckit::LocalConfiguration conf = createConfig();

  const eckit::LocalConfiguration sub(conf, "sub");
  EXPECT(util::isSubConfig(sub));

  const eckit::LocalConfiguration fp(conf, "fp");
  EXPECT_NOT(util::isSubConfig(fp));

  const eckit::LocalConfiguration vector(conf, "vector");
  EXPECT_NOT(util::isSubConfig(vector));
}

CASE("util/ConfigFunctions/isFinal") {
  eckit::LocalConfiguration conf = createConfig();

  const eckit::LocalConfiguration sub(conf, "sub");
  EXPECT_NOT(util::isFinal(sub));

  const eckit::LocalConfiguration fp(conf, "fp");
  EXPECT(util::isFinal(fp));

  const eckit::LocalConfiguration vector(conf, "vector");
  EXPECT_NOT(util::isFinal(vector));
}

CASE("util/ConfigFunctions/seekAndReplaceString") {
  eckit::LocalConfiguration conf = createConfig();

  eckit::LocalConfiguration sub(conf, "sub");
  util::seekAndReplace(sub, "%pattern%", "001");
  std::stringstream ss;
  ss << sub << std::endl;
  const std::string str = ss.str();

  const eckit::LocalConfiguration ref(conf, "ref");
  std::stringstream ssRef;
  ssRef << ref << std::endl;
  const std::string strRef = ssRef.str();

  EXPECT(str == strRef);
}

CASE("util/ConfigFunctions/seekAndReplaceIndex") {
  eckit::LocalConfiguration conf = createConfig();

  eckit::LocalConfiguration sub(conf, "sub");
  util::seekAndReplace(sub, "%pattern%", 1, 3);
  std::stringstream ss;
  ss << sub << std::endl;
  const std::string str = ss.str();

  const eckit::LocalConfiguration ref(conf, "ref");
  std::stringstream ssRef;
  ssRef << ref << std::endl;
  const std::string strRef = ssRef.str();

  EXPECT(str == strRef);
}

CASE("util/ConfigFunctions/mergeConfigs") {
  eckit::LocalConfiguration conf = createConfig();

  const eckit::LocalConfiguration ref(conf, "ref");
  const eckit::LocalConfiguration other(conf, "other");
  eckit::LocalConfiguration config = util::mergeConfigs(ref, other);
  std::stringstream ss;
  ss << config << std::endl;
  const std::string str = ss.str();

  const eckit::LocalConfiguration merged(conf, "merged");
  std::stringstream ssMerged;
  ssMerged << merged << std::endl;
  const std::string strMerged = ssMerged.str();

  EXPECT(str == strMerged);
}

class ConfigFunctions : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "test::ConfigFunctions";}

  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test

#endif  // TEST_UTIL_CONFIGFUNCTIONS_H_
