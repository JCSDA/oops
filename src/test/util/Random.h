/*
 * (C) Copyright 2017 UCAR
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 *
 */

#ifndef TEST_UTIL_RANDOM_H_
#define TEST_UTIL_RANDOM_H_

#define BOOST_TEST_NO_MAIN
#define BOOST_TEST_ALTERNATIVE_INIT_API
#define BOOST_TEST_DYN_LINK

#include <iostream>
#include <string>
#include <boost/test/unit_test.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/runs/Test.h"
#include "oops/util/Random.h"
#include "test/TestEnvironment.h"

namespace test {

// -----------------------------------------------------------------------------

void testCppRandom() {  

  const eckit::LocalConfiguration config(TestEnvironment::config(),"test_random");
  std::size_t N = static_cast<size_t>(config.getInt("N"));
  unsigned int seed = static_cast<unsigned int>(config.getInt("seed"));
  
  BOOST_CHECK_EQUAL(N, 10);  
  BOOST_CHECK_EQUAL(seed, 1234321);  

  std::vector<float> real_range = config.getFloatVector("uniform_real_range");
  util::UniformDistribution<float> x(N,real_range[0],real_range[1],seed);
  std::cout << "\nMSM Uniform: \n" << x << std::endl;
  
  std::vector<int> int_range = config.getIntVector("uniform_int_range");
  util::UniformIntDistribution<int> y(N,int_range[0],int_range[1],seed);
  std::cout << "\nMSM Int: \n" << y << std::endl;

  double normal_mean = config.getDouble("normal_mean");
  double normal_sdev = config.getDouble("normal_sdev");
  util::NormalDistribution<double> z(N,normal_mean,normal_sdev,seed);
  std::cout << "\nMSM Normal: \n" << z << std::endl;
  
}

// -----------------------------------------------------------------------------

class Random : public oops::Test {

 public:
  Random() {}
  virtual ~Random() {}

 private:
  std::string testid() const {return "test::Random";}

  void register_tests() const {
    boost::unit_test::test_suite * ts = BOOST_TEST_SUITE("util/Random");

    ts->add(BOOST_TEST_CASE(&testCppRandom));

    boost::unit_test::framework::master_test_suite().add(ts);
  }  
};

// -----------------------------------------------------------------------------
    
} // namespace test

#endif  // TEST_UTIL_RANDOM_H_

