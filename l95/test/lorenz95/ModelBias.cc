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
#include <boost/shared_ptr.hpp>

#include "./TestConfig.h"
#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelBiasCorrection.h"
#include "lorenz95/ModelBiasCovariance.h"
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
    biasconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "ModelBias"));
    nobias_.reset(new eckit::LocalConfiguration());
    eckit::LocalConfiguration covarCfg(TestConfig::config(), "Covariance");
    lorenz95::ModelBiasCovariance covar(covarCfg, *resol_);
    bias1_ = biasconf_->getDouble("bias");
    bias2_ = 2.5 * bias1_;
    dbias_.reset(new lorenz95::ModelBiasCorrection(*resol_, covar.config()));
    dbias_->bias() = bias2_;
  }
  ~ModBiasTestFixture() {}
  boost::scoped_ptr<lorenz95::Resolution> resol_;
  boost::scoped_ptr<const eckit::LocalConfiguration> biasconf_;
  boost::scoped_ptr<const eckit::LocalConfiguration> nobias_;
  boost::scoped_ptr<lorenz95::ModelBiasCorrection> dbias_;
  double bias1_;
  double bias2_;
};
// -----------------------------------------------------------------------------
CASE("test_modelBias") {
  ModBiasTestFixture f;
// -----------------------------------------------------------------------------
  SECTION("test_modelBias_constructor_bias") {
    lorenz95::ModelBias mbias(*f.resol_, *f.biasconf_);
    EXPECT(mbias.bias() == f.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBias_constructor_nobias") {
    lorenz95::ModelBias mbias(*f.resol_, *f.nobias_);

    // because the biasconf_ is empty when used,
    // the active_ flag is false and the bias_ value is 0.0
    EXPECT(mbias.bias() == 0.0);
}
// -----------------------------------------------------------------------------
  SECTION("test_modelBias_constructor_active") {
    lorenz95::ModelBias mbias1(*f.resol_, *f.biasconf_);
    lorenz95::ModelBias mbias2(*f.resol_, mbias1);
    EXPECT(mbias2.bias() == f.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBias_constructor_inactive") {
    lorenz95::ModelBias mbias1(*f.resol_, *f.nobias_);
    lorenz95::ModelBias mbias2(*f.resol_, mbias1);

    // because the biasconf_ is empty when used,
    // the active_ flag is false and the bias_ value is 0.0
    EXPECT(mbias2.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBias_constructor_active_copy") {
    lorenz95::ModelBias mbias1(*f.resol_, *f.biasconf_);
    lorenz95::ModelBias mbias2(mbias1, true);
    EXPECT(mbias2.bias() == f.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBias_constructor_active_nocopy") {
    lorenz95::ModelBias mbias1(*f.resol_, *f.biasconf_);
    lorenz95::ModelBias mbias2(mbias1, false);

    // because copy flag is set to false, bias is 0.0
    EXPECT(mbias2.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBias_constructor_inactive_copy") {
    lorenz95::ModelBias mbias1(*f.resol_, *f.nobias_);
    lorenz95::ModelBias mbias2(mbias1, true);

    // because the biasconf_ is empty when used,
    // the active_ flag is false and the bias_ value is 0.0
    EXPECT(mbias2.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBias_constructor_inactive_nocopy") {
    lorenz95::ModelBias mbias1(*f.resol_, *f.nobias_);
    lorenz95::ModelBias mbias2(mbias1, false);
    EXPECT(mbias2.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBias_destructor") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBias_classname") {
    lorenz95::ModelBias mbias(*f.resol_, *f.biasconf_);
    EXPECT(mbias.classname() == "lorenz95::ModelBias");
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBias_compound_assignment_active") {
    lorenz95::ModelBias mbias(*f.resol_, *f.biasconf_);

    mbias += *f.dbias_;

    EXPECT(mbias.bias() == f.bias1_ + f.bias2_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBias_compound_assignment_inactive") {
    lorenz95::ModelBias mbias(*f.resol_, *f.nobias_);

    mbias += *f.dbias_;

    // because the biasconf_ is empty,
    // mbias bias value is not modified from 0.0
    EXPECT(mbias.bias() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBias_bias") {
    lorenz95::ModelBias mbias(*f.resol_, *f.biasconf_);
    EXPECT(mbias.bias() == f.bias1_);
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBias_read") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBias_write") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_modelBias_stream_output") {
    lorenz95::ModelBias mbias(*f.resol_, *f.biasconf_);

    // use the operator<< method to write the value to a file
    std::filebuf fb;
    std::string filename("ModelBiasTest.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << mbias;
    fb.close();

    // then read the value that was written to the file
    std::string inputString;
    std::string inputBias;
    double testBias = f.bias1_;
    double bias = 0.0;
    int biasStartPos = 12;  // length of "ModelBias = " is 12
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
}  //  CASE
// -----------------------------------------------------------------------------
}  // namespace test
int main(int argc, char **argv)
{
    return eckit::testing::run_tests ( argc, argv );
}
