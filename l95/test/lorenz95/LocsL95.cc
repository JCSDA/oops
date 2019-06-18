/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <memory>  //  for std::unique_ptr


#include "./TestConfig.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "lorenz95/LocsL95.h"
#include "lorenz95/ObsTable.h"
#include "oops/runs/Test.h"
#include "oops/util/DateTime.h"
#include "test/TestFixture.h"

namespace test {

// -----------------------------------------------------------------------------
class LocsTestFixture : TestFixture {
 public:
  LocsTestFixture() {
    const eckit::LocalConfiguration conf(TestConfig::config(), "Observations");
    const util::DateTime bgn(conf.getString("window_begin"));
    const util::DateTime end(conf.getString("window_end"));
    const eckit::LocalConfiguration otconf(conf, "Observation.ObsSpace");
    ot_.reset(new lorenz95::ObsTable(otconf, bgn, end));
    t1_.reset(new util::DateTime("2010-01-01T12:00:00Z"));
    t2_.reset(new util::DateTime("2010-01-02T00:00:00Z"));
  }
  ~LocsTestFixture() {}
  std::unique_ptr<lorenz95::ObsTable> ot_;
  std::unique_ptr<util::DateTime> t1_;
  std::unique_ptr<util::DateTime> t2_;
};
// -----------------------------------------------------------------------------
CASE("test_LocsL95") {
  LocsTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_LocsL95_constructor") {
    std::unique_ptr<lorenz95::LocsL95> locs(fix.ot_->locations(*fix.t1_, *fix.t2_));
    EXPECT(locs.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_LocsL95_nobs") {
    std::unique_ptr<lorenz95::LocsL95> locs(fix.ot_->locations(*fix.t1_, *fix.t2_));
    size_t ref = 80;
    EXPECT(locs->size() == ref);
  }
// -----------------------------------------------------------------------------
  SECTION("test_LocsL95_operator") {
    std::unique_ptr<lorenz95::LocsL95> locs(fix.ot_->locations(*fix.t1_, *fix.t2_));
    double pos = 0.0;
    for (size_t jj = 0; jj < locs->size(); ++jj) {
      EXPECT(oops::is_close((*locs)[jj], pos, 0.00000001));

      pos += 0.05;
      if (pos >= 1.0) pos=0.0;
    }
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
