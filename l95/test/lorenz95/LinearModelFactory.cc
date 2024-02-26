/*
 * (C) Copyright 2020 Met Office UK.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "lorenz95/L95Traits.h"
#include "lorenz95/LocsL95.h"
#include "oops/generic/LinearModelBase.h"
#include "oops/runs/Run.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "test/TestEnvironment.h"

namespace test {

CASE("test_linearmodelfactory") {
  EXPECT_EQUAL(oops::LinearModelFactory<lorenz95::L95Traits>::getMakerNames(),
               std::vector<std::string>{"L95TLM"});
}

class LinearModelFactory : public oops::Test {
 public:
  LinearModelFactory() {
    // Constructing a LocsL95 object introduces a dependency on a non-inline function from
    // the lorenz95 library and thus prevents the linker from dropping that library from the list of
    // dependencies. This ensures static initialization of that library is completed before main()
    // starts and hence the L95TLM linear model is registered in the linear model factory.
    lorenz95::LocsL95 locs({}, {});
  }

 private:
  std::string testid() const override {return "test::LinearModelFactory";}

  void register_tests() const override {}

  void clear() const override {}
};

}  // namespace test

int main(int argc, char **argv) {
  oops::Run run(argc, argv);
  test::LinearModelFactory tests;
  return run.execute(tests);
}
