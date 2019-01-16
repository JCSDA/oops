/*
 * (C) Copyright 2019 UCAR
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

#include <cmath>
#include <iostream>
#include <limits>
#include <string>
#include <vector>
#include <boost/noncopyable.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/test/unit_test.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/runs/Test.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"
#include "test/TestEnvironment.h"
#include "test/util/Fortran.h"

namespace test {

// -----------------------------------------------------------------------------
/*! Test Fixture for random number generators */

class RandomFixture : private boost::noncopyable {
 public:
  static const eckit::Configuration & test() {return *getInstance().test_;}

 private:
  static RandomFixture & getInstance() {
    static RandomFixture theRandomFixture;
    return theRandomFixture;
  }

  RandomFixture() {
    test_.reset(new eckit::LocalConfiguration(TestEnvironment::config(), "test_random"));
  }

  ~RandomFixture() {}

  boost::scoped_ptr<const eckit::LocalConfiguration> test_;
};

// -----------------------------------------------------------------------------
/*! Test C++ implementation of random number generators
 *
 * \details This is intended to make sure that the random number generators in 
 * oops::util::Random.h are reproducible.  Specifically, we want to make sure that 
 * they are independent of the compiler and platform for a given random seed.
 *
*/

void testCppRandom() {
  typedef RandomFixture Test_;
  
  std::size_t N = static_cast<size_t>(Test_::test().getInt("N"));
  unsigned int seed = static_cast<unsigned int>(Test_::test().getInt("seed"));

  /*! Test uniform real distrubution 
   * The tolerance is based on the precision of the data type, in this case <float>.  
   * However, we have to multiply this by 100 because boost wants the tolerance as 
   * a percentage.  We also have to take into account the magnitude of the number 
   * itself, which requires another multiplication.  This can potentially decrease 
   * the accuracy by one significant digit.  So, also include a safety factor sfac.
  */
  double sfac = 10;
  float tol = 100 * sfac * std::numeric_limits<float>::epsilon();
  std::vector<float> real_range = Test_::test().getFloatVector("uniform_real_range");
  util::UniformDistribution<float> x(N, real_range[0], real_range[1], seed);
  std::vector<float> x_check = Test_::test().getFloatVector("uniform_real_answer");
  oops::Log::info() << "\nTesting oops::util::Random.h Uniform Real Distribution: \n"
                    << x << std::endl;
  for (std::size_t jj = 0; jj < N; ++jj) BOOST_CHECK_CLOSE(x[jj], x_check[jj],
                                                           tol * std::abs(x_check[jj]));

  /*! Test uniform integer distribution */
  std::vector<int> int_range = Test_::test().getIntVector("uniform_int_range");
  util::UniformIntDistribution<int> y(N, int_range[0], int_range[1], seed);
  std::vector<int> y_check = Test_::test().getIntVector("uniform_int_answer");
  oops::Log::info() << "\nTesting oops::util::Random.h Uniform Int Distribution: \n"
                    << y << std::endl;
  for (std::size_t jj = 0; jj < N; ++jj) BOOST_CHECK_EQUAL(y[jj], y_check[jj]);

  /*! Test normal distribution */
  tol = 100 * sfac * std::numeric_limits<double>::epsilon();  // see comment above
  double normal_mean = Test_::test().getDouble("normal_mean");
  double normal_sdev = Test_::test().getDouble("normal_sdev");
  util::NormalDistribution<double> z(N, normal_mean, normal_sdev, seed);
  std::vector<double> z_check = Test_::test().getDoubleVector("normal_answer");
  oops::Log::info() << "\nTesting oops::util::Random.h Gaussian Distribution: \n"
                    << z << std::endl;
  for (std::size_t jj = 0; jj < N; ++jj) BOOST_CHECK_CLOSE(z[jj], z_check[jj],
                                                           tol * std::abs(z_check[jj]));
  }

// -----------------------------------------------------------------------------
/*! Test Fortran implementation of random number generators
 *
 * \details This is intended to make sure that the random number generators in 
 * oops::util::Random.h are reproducible.  Specifically, we want to make sure that 
 * they are independent of the compiler and platform for a given random seed.
 *
 */

void testFortranRandom() {
  typedef RandomFixture Test_;

  const eckit::Configuration * config = &Test_::test();

  BOOST_CHECK(test_uniform_real_f(&config) == 0);
  
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
    ts->add(BOOST_TEST_CASE(&testFortranRandom));

    boost::unit_test::framework::master_test_suite().add(ts);
  }
};

// -----------------------------------------------------------------------------
}  // namespace test

#endif  // TEST_UTIL_RANDOM_H_

