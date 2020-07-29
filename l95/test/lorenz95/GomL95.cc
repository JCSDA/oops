/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <memory>

#include "./TestConfig.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "lorenz95/GomL95.h"
#include "lorenz95/LocsL95.h"
#include "lorenz95/ObservationL95.h"
#include "lorenz95/ObsTable.h"
#include "oops/base/Variables.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class GomTestFixture : TestFixture {
 public:
  GomTestFixture() {
    const eckit::LocalConfiguration conf(TestConfig::config(), "observations");
    const util::DateTime bgn(conf.getString("window begin"));
    const util::DateTime end(conf.getString("window end"));
    const eckit::LocalConfiguration otconf(conf, "observation.obs space");
    lorenz95::ObsTable ot(otconf, oops::mpi::comm(), bgn, end);
    util::DateTime t1("2010-01-01T03:00:00Z");
    util::DateTime t2("2010-01-02T06:00:00Z");
    locs_.reset(ot.locations(t1, t2));
    novar_.reset(new oops::Variables());
  }
  ~GomTestFixture() {}
  std::unique_ptr<lorenz95::LocsL95> locs_;
  std::unique_ptr<oops::Variables> novar_;
};
// -----------------------------------------------------------------------------
CASE("test_GomL95") {
  GomTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_GomL95_constructor") {
    std::unique_ptr<lorenz95::GomL95> gom(new lorenz95::GomL95(*fix.locs_, *fix.novar_));
    EXPECT(gom.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_GomL95_nobs") {
    std::unique_ptr<lorenz95::GomL95> gom(new lorenz95::GomL95(*fix.locs_, *fix.novar_));
    size_t ref = 160;
    EXPECT(gom->size() == ref);
  }
// -----------------------------------------------------------------------------
  SECTION("test_gomL95_classname") {
    std::unique_ptr<lorenz95::GomL95> gom(new lorenz95::GomL95(*fix.locs_, *fix.novar_));
    EXPECT(gom->classname() == "lorenz95::GomL95");
  }
// -----------------------------------------------------------------------------
  SECTION("test_gomL95_zero") {
    std::unique_ptr<lorenz95::GomL95> gom(new lorenz95::GomL95(*fix.locs_, *fix.novar_));
    gom->zero();
    for (size_t i = 0; i < gom->size(); ++i) {
      EXPECT((*gom)[i] == 0.0);
    }
  }
// -----------------------------------------------------------------------------
  SECTION("test_gomL95_dot_product_with") {
    std::unique_ptr<lorenz95::GomL95> gom1(new lorenz95::GomL95(*fix.locs_, *fix.novar_));
    gom1->zero();
    std::unique_ptr<lorenz95::GomL95> gom2(new lorenz95::GomL95(*fix.locs_, *fix.novar_));
    gom2->zero();

    double zz = gom1->dot_product_with(*gom2);
    EXPECT(zz == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_gomL95_operator") {
    std::unique_ptr<lorenz95::GomL95> gom1(new lorenz95::GomL95(*fix.locs_, *fix.novar_));
    gom1->zero();
    std::unique_ptr<lorenz95::GomL95> gom2(new lorenz95::GomL95(*fix.locs_, *fix.novar_));
    gom2->zero();

    (*gom1)[1] = 1.0;
    (*gom2)[3] = 1.0;
    double zz = gom1->dot_product_with(*gom2);
    EXPECT(zz == 0.0);

    (*gom1)[3] = 1.0;
    zz = gom1->dot_product_with(*gom2);
    EXPECT(zz == 1.0);
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------

}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
