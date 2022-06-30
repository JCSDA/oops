/*
 * (C) Crown copyright 2021 Met Office UK.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/generic/ObsErrorDiag.h"
#include "oops/runs/Run.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "test/TestEnvironment.h"

namespace test {

CASE("test_obserrordiagparameters") {
  for (const eckit::LocalConfiguration &conf :
       TestEnvironment::config().getSubConfigurations("cases")) {
    eckit::LocalConfiguration obsErrorConf(conf, "obs error");
    const std::string expectedMsg = conf.getString("expected error");
    oops::ObsErrorDiagParameters parameters;
    EXPECT_THROWS_MSG(parameters.deserialize(obsErrorConf), expectedMsg.c_str());
  }
}

class ObsErrorDiagParameters : public oops::Test {
 public:
  ObsErrorDiagParameters() {}

 private:
  std::string testid() const override {return "test::ObsErrorDiagParameters";}

  void register_tests() const override {}

  void clear() const override {}
};

}  // namespace test

int main(int argc, char **argv) {
  oops::Run run(argc, argv);
  test::ObsErrorDiagParameters tests;
  return run.execute(tests);
}
