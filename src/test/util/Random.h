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

#include "eckit/config/LocalConfiguration.h"
#include <boost/test/unit_test.hpp>
#include "oops/runs/Test.h"
#include "oops/util/Random.h"

namespace test {

// -----------------------------------------------------------------------------

void testCppRandom() {  

  const eckit::LocalConfiguration config(TestEnvironment::config(),"test_random");
  size_t N = static_cast<size_t>(config.getInt("N"));
  unsigned int seed = static_cast<unsigned int>(config.getInt("seed"));
  
  BOOST_CHECK_EQUAL(seed, 1234321);  
  
  //util::UniformDistribution<int> x(N,1,100,seed);

  //for {size_t jj=0; jj < N; ++j) {
  //  std::cout << "MSM Random: " << jj << ": ", x[jj] << std::endl;
  //}
  
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

