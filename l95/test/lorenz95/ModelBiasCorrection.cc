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
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelBiasCorrection.h"
#include "lorenz95/Resolution.h"
#include "test/TestFixture.h"

using eckit::types::is_approximately_equal;

namespace test {

// -----------------------------------------------------------------------------
class ModBiasTestFixture : TestFixtureBase<false> {
 public:
  ModBiasTestFixture() {
    eckit::LocalConfiguration res(TestConfig::config(), "resolution");
    resol_.reset(new lorenz95::Resolution(res));
    conf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "ModelBiasCovariance"));
    nobias_.reset(new eckit::LocalConfiguration());
    bias1_ = TestConfig::config().getDouble("ModelBias.bias");
    bias2_ = 2.5 * bias1_;
    fact_ = 1.2345;
  }
  ~ModBiasTestFixture() {}
  boost::scoped_ptr<lorenz95::Resolution> resol_;
  boost::scoped_ptr<const eckit::LocalConfiguration> conf_;
  boost::scoped_ptr<const eckit::LocalConfiguration> nobias_;
  double bias1_;
  double bias2_;
  double fact_;
};
// -----------------------------------------------------------------------------
CASE("test_modelBiasCorrection") {
  ModBiasTestFixture f;
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_constructor_config") {
    boost::scoped_ptr<lorenz95::ModelBiasCorrection> dx(
      new lorenz95::ModelBiasCorrection(*f.resol_, *f.conf_));
    EXPECT(dx.get() != NULL);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_constructor_no_config") {
    boost::scoped_ptr<lorenz95::ModelBiasCorrection> dx(
      new lorenz95::ModelBiasCorrection(*f.resol_, *f.nobias_));
    EXPECT(dx.get() != NULL);
}
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_copy_ctor_active_copy") {
    lorenz95::ModelBiasCorrection dx1(*f.resol_, *f.conf_);
    dx1.bias() = f.bias1_;

    lorenz95::ModelBiasCorrection dx2(dx1, true);

    EXPECT(dx2.bias() == f.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_copy_ctor_active_no_copy") {
    lorenz95::ModelBiasCorrection dx1(*f.resol_, *f.conf_);
    dx1.bias() = f.bias1_;

    // construct a copy of it with the copy flag set to false
    lorenz95::ModelBiasCorrection dx2(dx1, false);

    // because the copy is false,
    // the active_ flag is true and the bias_ value is 0.0
    EXPECT(dx2.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_copy_ctor_inactive_copy") {
    lorenz95::ModelBiasCorrection dx1(*f.resol_, *f.nobias_);
    dx1.bias() = f.bias1_;

    // construct a copy of it with the copy flag set to true
    lorenz95::ModelBiasCorrection dx2(dx1, true);

    // because the cfg is empty when used,
    // the active_ flag is false and the bias_ value is 0.0
    EXPECT(dx2.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_copy_ctor_inactive_no_copy") {
    lorenz95::ModelBiasCorrection dx1(*f.resol_, *f.nobias_);
    dx1.bias() = f.bias1_;

    // construct a copy of it with the copy flag set to false
    lorenz95::ModelBiasCorrection dx2(dx1, false);

    // because the cfg is empty when used and the copy flag is false,
    // the active_ flag is false and the bias_ value is 0.0
    EXPECT(dx2.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_copy_ctor_config_active") {
    lorenz95::ModelBiasCorrection dx1(*f.resol_, *f.conf_);
    dx1.bias() = f.bias1_;

    lorenz95::ModelBiasCorrection dx2(dx1, *f.conf_);

    EXPECT(dx2.bias() == f.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_copy_ctor_config_inactive") {
    lorenz95::ModelBiasCorrection dx1(*f.resol_, *f.nobias_);
    dx1.bias() = f.bias1_;

    lorenz95::ModelBiasCorrection dx2(dx1, *f.conf_);

     // because the covarCfg is empty when used (regardless of the cfg),
     // the active_ flag is false and the bias_ value is 0.0
     EXPECT(dx2.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_diff_active") {
    lorenz95::ModelBias xx1(*f.resol_, *f.conf_);
    xx1.bias() = f.bias1_;
    lorenz95::ModelBias xx2(*f.resol_, *f.conf_);
    xx2.bias() = f.bias2_;

    lorenz95::ModelBiasCorrection dx(*f.resol_, *f.conf_);

    dx.diff(xx1, xx2);

    EXPECT(dx.bias() == f.bias1_ - f.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_diff_inactive") {
    lorenz95::ModelBias xx1(*f.resol_, *f.conf_);
    xx1.bias() = f.bias1_;
    lorenz95::ModelBias xx2(*f.resol_, *f.conf_);
    xx2.bias() = f.bias2_;

    lorenz95::ModelBiasCorrection dx(*f.resol_, *f.nobias_);

    dx.diff(xx1, xx2);

    // because the active_ flag is false, the bias cannot be updated
    EXPECT(dx.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_zero") {
    lorenz95::ModelBiasCorrection dx(*f.resol_, *f.conf_);
    dx.bias() = f.bias1_;

    dx.zero();

    EXPECT(dx.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_assignment_active") {
    lorenz95::ModelBiasCorrection dx1(*f.resol_, *f.conf_);
    dx1.bias() = f.bias1_;
    lorenz95::ModelBiasCorrection dx2(*f.resol_, *f.conf_);
    dx2.bias() = f.bias2_;

    dx1 = dx2;

    // the original MBC should have the same bias value as the copy MBC
    EXPECT(dx1.bias() == f.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_assignment_inactive") {
    lorenz95::ModelBiasCorrection dx1(*f.resol_, *f.nobias_);
    dx1.bias() = f.bias1_;

    lorenz95::ModelBiasCorrection dx2(*f.resol_, *f.conf_);
    dx2.bias() = f.bias2_;

    dx1 = dx2;

    // the active_ value is zero, so the bias will be zero
    EXPECT(dx1.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_assignment_add_active") {
    lorenz95::ModelBiasCorrection dx1(*f.resol_, *f.conf_);
    dx1.bias() = f.bias1_;
    lorenz95::ModelBiasCorrection dx2(*f.resol_, *f.conf_);
    dx2.bias() = f.bias2_;

    dx1 += dx2;

    EXPECT(dx1.bias() == f.bias1_ + f.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_assignment_add_inactive") {
    lorenz95::ModelBiasCorrection dx1(*f.resol_, *f.nobias_);
    dx1.bias() = f.bias1_;
    lorenz95::ModelBiasCorrection dx2(*f.resol_, *f.nobias_);
    dx2.bias() = f.bias2_;

    dx1 += dx2;

    // the active_ value is zero, so the bias will be unchanged
    EXPECT(dx1.bias() == f.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_assignment_subtract_active") {
    lorenz95::ModelBiasCorrection dx1(*f.resol_, *f.conf_);
    dx1.bias() = f.bias1_;
    lorenz95::ModelBiasCorrection dx2(*f.resol_, *f.conf_);
    dx2.bias() = f.bias2_;

    dx1 -= dx2;

    EXPECT(dx1.bias() == f.bias1_ - f.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_assignment_subtract_inactive") {
    lorenz95::ModelBiasCorrection dx1(*f.resol_, *f.nobias_);
    dx1.bias() = f.bias1_;
    lorenz95::ModelBiasCorrection dx2(*f.resol_, *f.nobias_);
    dx2.bias() = f.bias2_;

    dx1 -= dx2;

    // the active_ value is zero, so the bias will be unchanged
    EXPECT(dx1.bias() == f.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_assignment_multiply_active") {
    lorenz95::ModelBiasCorrection dx(*f.resol_, *f.conf_);
    dx.bias() = f.bias1_;

    dx *= f.fact_;

    EXPECT(dx.bias() == f.bias1_ * f.fact_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_assignment_multiply_inactive") {
    lorenz95::ModelBiasCorrection dx(*f.resol_, *f.nobias_);
    dx.bias() = f.bias1_;

    dx *= f.fact_;

    // the active_ value is zero, so the bias will be unchanged
    EXPECT(dx.bias() == f.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_axpy_active") {
    lorenz95::ModelBiasCorrection dx1(*f.resol_, *f.conf_);
    dx1.bias() = f.bias1_;
    lorenz95::ModelBiasCorrection dx2(*f.resol_, *f.conf_);
    dx2.bias() = f.bias2_;

    dx1.axpy(f.fact_, dx2);

    EXPECT(dx1.bias() == (f.bias1_ + f.fact_ * f.bias2_));
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_axpy_inactive") {
    lorenz95::ModelBiasCorrection dx1(*f.resol_, *f.nobias_);
    dx1.bias() = f.bias1_;
    lorenz95::ModelBiasCorrection dx2(*f.resol_, *f.nobias_);
    dx2.bias() = f.bias2_;

    dx1.axpy(f.fact_, dx2);

    // the active_ value is zero, so the bias will be unchanged
    EXPECT(dx1.bias() == f.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_dot_product_with_active") {
    lorenz95::ModelBiasCorrection dx1(*f.resol_, *f.conf_);
    dx1.bias() = f.bias1_;
    lorenz95::ModelBiasCorrection dx2(*f.resol_, *f.conf_);
    dx2.bias() = f.bias2_;

    double dpwResult = dx1.dot_product_with(dx2);

    EXPECT(dpwResult == f.bias1_ * f.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_dot_product_with_inactive") {
    lorenz95::ModelBiasCorrection dx1(*f.resol_, *f.nobias_);
    dx1.bias() = f.bias1_;
    lorenz95::ModelBiasCorrection dx2(*f.resol_, *f.nobias_);
    dx2.bias() = f.bias2_;

    double dpwResult = dx1.dot_product_with(dx2);

    // because of the empty config, the result is 0.0
    EXPECT(dpwResult == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_read") {
    // nothing to test
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_write") {
    // nothing to test
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_stream_output") {
    lorenz95::ModelBiasCorrection dx(*f.resol_, *f.conf_);
    dx.bias() = f.bias1_;

    // use the operator<< method to write the value to a file
    std::filebuf fb;
    std::string filename("ModelBiasCorrectionTest.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << dx;
    fb.close();

    // then read the value that was written to the file
    std::string inputString;
    std::string inputBias;
    double testBias = f.bias1_;
    double bias = 0.0;
    int biasStartPos = 22;  // length of "ModelBiasCorrection = " is 22
    std::ifstream inputFile(filename.c_str());
    if (inputFile.is_open()) {
      getline(inputFile, inputString);  // ignore first (blank) line
      getline(inputFile, inputString);

      inputBias = inputString.substr(biasStartPos);

      try {
        bias = boost::lexical_cast<double>(inputBias);
      }
      catch(boost::bad_lexical_cast const&) {
        std::cout << "operator<< incorrectly output a non-double" << std::endl;
      }

      EXPECT(is_approximately_equal(testBias, bias, 0.0001));
    } else {
      // if we can't open the file then we can't
      // verify that the value was correctly written
      std::cout << "operator<< functionality cannot be determined" << std::endl;
    }
    inputFile.close();
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBiasCorrection_bias") {
    lorenz95::ModelBiasCorrection dx(*f.resol_, *f.conf_);
    dx.bias() = f.bias1_;

    // this one test checks both the setting and getting of the bias value
    EXPECT(dx.bias() == f.bias1_);
  }
// -----------------------------------------------------------------------------
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
