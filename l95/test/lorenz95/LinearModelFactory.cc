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
#include "lorenz95/TLML95.h"
#include "oops/interface/LinearModel.h"
#include "oops/runs/Run.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "test/TestEnvironment.h"

namespace test {

CASE("test_linearmodelparameterswrapper_valid_name") {
  eckit::LocalConfiguration config(TestEnvironment::config(), "valid linear model name");
  oops::LinearModelParametersWrapper<lorenz95::L95Traits> parameters;
  EXPECT_NO_THROW(parameters.validateAndDeserialize(config));

  // Constructing a TLML95 object here introduces a dependency on a non-inline function from
  // the lorenz95 library and thus prevents the linker from dropping that library from the list of
  // dependencies. This ensures static initialization of that library is completed before main()
  // starts and hence the L95TLM linear model is registered in the linear model factory.
  lorenz95::Resolution resol(40);
  lorenz95::TLML95 tlm(resol, config);
}

CASE("test_linearmodelparameterswrapper_invalid_name") {
  eckit::LocalConfiguration config(TestEnvironment::config(), "invalid linear model name");
  oops::LinearModelParametersWrapper<lorenz95::L95Traits> parameters;
  if (oops::Parameters::isValidationSupported())
    EXPECT_THROWS_MSG(parameters.validate(config), "unrecognized enum value");
  EXPECT_THROWS_MSG(parameters.deserialize(config),
                    "does not exist in the tangent linear model factory");
}

CASE("test_linearmodelfactory") {
  EXPECT_EQUAL(oops::LinearModelFactory<lorenz95::L95Traits>::getMakerNames(),
               std::vector<std::string>{"L95TLM"});
}

class LinearModelFactory : public oops::Test {
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
