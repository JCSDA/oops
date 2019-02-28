/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <fstream>
#include <iostream>

#include <boost/lexical_cast.hpp>
#include <boost/scoped_ptr.hpp>

#include "./TestConfig.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "lorenz95/ObsBias.h"
#include "lorenz95/ObsBiasCorrection.h"
#include "oops/util/Logger.h"
#include "test/TestFixture.h"

using eckit::types::is_approximately_equal;

namespace test {

// -----------------------------------------------------------------------------
class ObsBiasTestFixture : TestFixture {
 public:
  ObsBiasTestFixture() {
    off_.reset(new eckit::LocalConfiguration());
    conf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "ObsBiasCovariance"));
    eckit::LocalConfiguration bconf(TestConfig::config(), "ObsBias");
    obias_.reset(new lorenz95::ObsBias(bconf));
    bias1_ = bconf.getDouble("bias");
    bias2_ = 3.5 * bias1_;
    fact_ = 1.234;
  }
  ~ObsBiasTestFixture() {}
  boost::scoped_ptr<const eckit::LocalConfiguration> off_;
  boost::scoped_ptr<const eckit::LocalConfiguration> conf_;
  double bias1_;
  double bias2_;
  double fact_;
  boost::scoped_ptr<lorenz95::ObsBias> obias_;
};

// -----------------------------------------------------------------------------
CASE("test_obsBiasCorrection") {
  ObsBiasTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_constructor_active") {
    lorenz95::ObsBiasCorrection dob(*fix.conf_);
    EXPECT(dob.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_constructor_inactive") {
    lorenz95::ObsBiasCorrection dob(*fix.off_);
    EXPECT(dob.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_copy_ctor_active_copy") {
    lorenz95::ObsBiasCorrection dob1(*fix.conf_);
    dob1.value() = fix.bias1_;

    lorenz95::ObsBiasCorrection dob2(dob1, true);

    EXPECT(dob2.value() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_copy_ctor_active_no_copy") {
    lorenz95::ObsBiasCorrection dob1(*fix.conf_);
    dob1.value() = fix.bias1_;

    lorenz95::ObsBiasCorrection dob2(dob1, false);

    // because the copy is false,
    // the active_ flag is true and the bias1_ value is 0.0
    EXPECT(dob2.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_copy_ctor_inactive_copy") {
    lorenz95::ObsBiasCorrection dob1(*fix.off_);
    dob1.value() = fix.bias1_;

    lorenz95::ObsBiasCorrection dob2(dob1, true);

    // because the cfg is empty when used,
    // the active_ flag is false and the bias1_ value is 0.0
    EXPECT(dob2.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_copy_ctor_inactive_no_copy") {
    lorenz95::ObsBiasCorrection dob1(*fix.off_);
    dob1.value() = fix.bias1_;

    lorenz95::ObsBiasCorrection dob2(dob1, false);

    // because the cfg is empty when used and the copy flag is false,
    // the active_ flag is false and the bias1_ value is 0.0
    EXPECT(dob2.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_copy_constructor_config") {
    lorenz95::ObsBiasCorrection dob1(*fix.conf_);
    dob1.value() = fix.bias1_;

    lorenz95::ObsBiasCorrection dob2(dob1, *fix.conf_);

    EXPECT(dob2.value() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_copy_constructor_no_config") {
    lorenz95::ObsBiasCorrection dob1(*fix.off_);
    dob1.value() = fix.bias1_;

    lorenz95::ObsBiasCorrection dob2(dob1, eckit::LocalConfiguration());

    // because the covarCfg is empty when used (regardless of the cfg),
    // the active_ flag is false and the bias1_ value is 0.0
    EXPECT(dob2.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_diff_active") {
    lorenz95::ObsBias obias2(*fix.conf_);
    obias2.value() = fix.bias2_;

    lorenz95::ObsBiasCorrection dob(*fix.conf_);

    dob.diff(*fix.obias_, obias2);

    EXPECT(dob.value() == fix.bias1_ - fix.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_diff_inactive") {
    lorenz95::ObsBias obias2(*fix.conf_);
    obias2.value() = fix.bias2_;

    // construct the dob object with empty config
    lorenz95::ObsBiasCorrection dob(*fix.off_);

    dob.diff(*fix.obias_, obias2);

    // because the OBC has empty config the active_flag is false and
    // the diff will not be performed, leaving the OBC value unchanged
    EXPECT(dob.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_zero") {
    lorenz95::ObsBiasCorrection dob(*fix.conf_);
    dob.value() = fix.bias1_;

    dob.zero();

    EXPECT(dob.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_assignment_active") {
    lorenz95::ObsBiasCorrection dob1(*fix.conf_);
    dob1.value() = fix.bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = fix.bias2_;

    dob1 = dob2;

    EXPECT(dob1.value() == fix.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_assignment_inactive") {
    lorenz95::ObsBiasCorrection dob1(*fix.off_);
    dob1.value() = fix.bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = fix.bias2_;

    dob1 = dob2;

    // because the OBC has empty config the active_ flag is false
    EXPECT(dob1.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_assignment_add_active") {
    lorenz95::ObsBiasCorrection dob1(*fix.conf_);
    dob1.value() = fix.bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = fix.bias2_;

    dob1 += dob2;

    EXPECT(dob1.value() == fix.bias1_ + fix.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_assign_add_inactive") {
    lorenz95::ObsBiasCorrection dob1(*fix.off_);
    dob1.value() = fix.bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = fix.bias2_;

    dob1 += dob2;

    // because the OBC has empty config, the bias value is unchanged
    EXPECT(dob1.value() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_assignment_subtract_active") {
    lorenz95::ObsBiasCorrection dob1(*fix.conf_);
    dob1.value() = fix.bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = fix.bias2_;

    dob1 -= dob2;

    EXPECT(dob1.value() == fix.bias1_ - fix.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_assignment_subtract_inactive") {
    lorenz95::ObsBiasCorrection dob1(*fix.off_);
    dob1.value() = fix.bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = fix.bias2_;

    dob1 -= dob2;

    // because the OBC has empty config, the bias value is unchanged
    EXPECT(dob1.value() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_assignment_multiply_active") {
    lorenz95::ObsBiasCorrection dob(*fix.conf_);
    dob.value() = fix.bias1_;

    dob *= fix.fact_;

    EXPECT(dob.value() == fix.bias1_ * fix.fact_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_assignment_multiply_inactive") {
    lorenz95::ObsBiasCorrection dob(*fix.off_);
    dob.value() = fix.bias1_;

    dob *= fix.fact_;

    // because the OBC has empty config, the bias value is unchanged
    EXPECT(dob.value() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_axpy_active") {
    lorenz95::ObsBiasCorrection dob1(*fix.conf_);
    dob1.value() = fix.bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = fix.bias2_;

    dob1.axpy(fix.fact_, dob2);

    EXPECT(dob1.value() == fix.bias1_ + (fix.fact_ * fix.bias2_));
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_axpy_inactive") {
    lorenz95::ObsBiasCorrection dob1(*fix.off_);
    dob1.value() = fix.bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = fix.bias2_;

    dob1.axpy(fix.fact_, dob2);

    // because the OBC has empty config, the bias value is unchanged
    EXPECT(dob1.value() == fix.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_dot_product_with_active") {
    lorenz95::ObsBiasCorrection dob1(*fix.conf_);
    dob1.value() = fix.bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = fix.bias2_;

    double dpwResult = dob1.dot_product_with(dob2);

    EXPECT(dpwResult == fix.bias1_ * fix.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_dot_product_with_inactive") {
    lorenz95::ObsBiasCorrection dob1(*fix.off_);
    dob1.value() = fix.bias1_;
    lorenz95::ObsBiasCorrection dob2(dob1);
    dob2.value() = fix.bias2_;

    double dpwResult = dob1.dot_product_with(dob2);

    // because the OBC has empty config, the result is 0
    EXPECT(dpwResult == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_read") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_write") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBiasCorrection_stream_output") {
    lorenz95::ObsBiasCorrection dob(*fix.conf_);
    dob.value() = fix.bias1_;

    // use the operator<< method to write the value to a file
    std::filebuf fb;
    std::string filename("ObsBiasCorrectionTest.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << dob;
    fb.close();

    // then read the value that was written to the file
    std::string inputString;
    std::string inputBias;
    double bias = 0.0;
    int biasStartPos = 20;  // length of "ObsBiasCorrection = " is 20
    std::ifstream inputFile(filename.c_str());
    if (inputFile.is_open()) {
      getline(inputFile, inputString);  // ignore first (blank) line
      getline(inputFile, inputString);

      inputBias = inputString.substr(biasStartPos);

      try {
        bias = boost::lexical_cast<double>(inputBias);
      }
      catch(boost::bad_lexical_cast const&) {
        oops::Log::error() << "operator<< incorrectly output a non-double" << std::endl;
      }

      EXPECT(is_approximately_equal(fix.bias1_, bias, 0.0001));
    } else {
      // if we can't open the file then we can't
      // verify that the value was correctly written
      oops::Log::error() << "operator<< functionality cannot be determined" << std::endl;
    }
    inputFile.close();
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
