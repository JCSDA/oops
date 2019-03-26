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
    biasconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "ObsBias"));
    nobias_.reset(new eckit::LocalConfiguration());
    covconf_.reset(new eckit::LocalConfiguration(TestConfig::config(), "Covariance"));
  }
  ~ObsBiasTestFixture() {}
  boost::scoped_ptr<const eckit::LocalConfiguration> biasconf_;
  boost::scoped_ptr<const eckit::LocalConfiguration> nobias_;
  boost::scoped_ptr<const eckit::LocalConfiguration> covconf_;
};
// -----------------------------------------------------------------------------
CASE("test_obsBias") {
  ObsBiasTestFixture fix;
// -----------------------------------------------------------------------------
  SECTION("test_obsBias_constructor_bias") {
    lorenz95::ObsBias ob(*fix.biasconf_);
    EXPECT(ob.value() == fix.biasconf_->getDouble("bias"));
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBias_constructor_no_bias") {
    lorenz95::ObsBias ob(*fix.nobias_);
    EXPECT(ob.value() == 0.0);
}
// -----------------------------------------------------------------------------
  SECTION("test_obsBias_copy_constructor_config_copy") {
    lorenz95::ObsBias ob1(*fix.biasconf_);
    lorenz95::ObsBias ob2(ob1, true);
    EXPECT(ob2.value() == fix.biasconf_->getDouble("bias"));
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBias_copy_constructor_config_no_copy") {
    lorenz95::ObsBias ob1(*fix.biasconf_);
    lorenz95::ObsBias ob2(ob1, false);

    // bias value is 0 because copy flag was set to false
    EXPECT(ob2.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBias_copy_constructor_no_config_copy") {
    lorenz95::ObsBias ob1(*fix.nobias_);
    lorenz95::ObsBias ob2(ob1, true);

    // bias value is 0 because an empty config was used
    EXPECT(ob2.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBias_copy_constructor_no_config_no_copy") {
    lorenz95::ObsBias ob1(*fix.nobias_);
    lorenz95::ObsBias ob2(ob1, true);

    // bias value is 0 because an empty config was used
    EXPECT(ob2.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBias_destructor") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBias_compound_assignment_add_active") {
    lorenz95::ObsBias ob(*fix.biasconf_);

    // construct an obsBiasCorrection object
    lorenz95::ObsBiasCorrection obsBiasCorrection(*fix.covconf_);
    obsBiasCorrection.value() = 3.14;

    ob += obsBiasCorrection;

    EXPECT(ob.value() == fix.biasconf_->getDouble("bias") + 3.14);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBias_compound_assignment_add_inactive") {
    lorenz95::ObsBias ob(*fix.nobias_);

    // construct an obsBiasCorrection object
    lorenz95::ObsBiasCorrection obsBiasCorrection(*fix.covconf_);
    obsBiasCorrection.value() = 3.14;

    ob += obsBiasCorrection;

    // the bias value will be 0 because the ob had an empty config
    EXPECT(ob.value() == 0.0);
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBias_value") {
    lorenz95::ObsBias ob(*fix.biasconf_);
    EXPECT(ob.value() == fix.biasconf_->getDouble("bias"));
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBias_read") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBias_write") {
    // not yet implemented
  }
// -----------------------------------------------------------------------------
  SECTION("test_obsBias_stream_output") {
    lorenz95::ObsBias ob(*fix.biasconf_);

    // use the operator<< method to write the value to a file
    std::filebuf fb;
    std::string filename("ObsBiasTest.txt");
    fb.open(filename.c_str(), std::ios::out);
    std::ostream os(&fb);
    os << ob;
    fb.close();

    // then read the value that was written to the file
    std::string inputString;
    std::string inputBias;
    double testBias = fix.biasconf_->getDouble("bias");
    double bias = 0.0;
    int biasStartPos = 10;  // length of "ObsBias = " is 10
    std::ifstream inputFile(filename.c_str());
    if (inputFile.is_open()) {
      getline(inputFile, inputString);
      getline(inputFile, inputString);

      inputBias = inputString.substr(biasStartPos);

      try {
        bias = boost::lexical_cast<double>(inputBias);
      }
      catch(boost::bad_lexical_cast const&) {
        oops::Log::error() << "operator<< incorrectly output a non-double" << std::endl;
      }

      EXPECT(is_approximately_equal(testBias, bias, 0.0001));
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
