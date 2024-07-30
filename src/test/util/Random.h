/*
 * (C) Copyright 2019 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#ifndef TEST_UTIL_RANDOM_H_
#define TEST_UTIL_RANDOM_H_

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <vector>

#define ECKIT_TESTING_SELF_REGISTER_CASES 0

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
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

  std::unique_ptr<const eckit::LocalConfiguration> test_;
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
   *  The tolerance is based on the precision of the data type, in this case <double>.
   *  We have to take into account the magnitude of the number
   *  itself, which requires a multiplication.  This can potentially decrease
   *  the accuracy by one significant digit.  So, also include a safety factor sfac.
  */
  double sfac = 10;
  double tol = sfac * std::numeric_limits<double>::epsilon();
  std::vector<double> real_range = Test_::test().getDoubleVector("uniform_real_range");
  util::UniformDistribution<double> x(N, real_range[0], real_range[1], seed);
  std::vector<double> x_check = Test_::test().getDoubleVector("uniform_real_answer");
  oops::Log::info() << "\nTesting oops::util::Random.h Uniform Real Distribution: \n"
                    << x << std::endl;
  for (std::size_t jj = 0; jj < N; ++jj)
    EXPECT(oops::is_close(x[jj], x_check[jj],
                                                 tol * std::abs(x_check[jj])));

  /*! Test uniform integer distribution */
  std::vector<int> int_range = Test_::test().getIntVector("uniform_int_range");
  util::UniformIntDistribution<int> y(N, int_range[0], int_range[1], seed);
  std::vector<int> y_check = Test_::test().getIntVector("uniform_int_answer");
  oops::Log::info() << "\nTesting oops::util::Random.h Uniform Int Distribution: \n"
                    << y << std::endl;
  for (std::size_t jj = 0; jj < N; ++jj) EXPECT(y[jj] == y_check[jj]);

  /*! Test normal distribution */
  tol = sfac * std::numeric_limits<double>::epsilon();  // see comment above
  double normal_mean = Test_::test().getDouble("normal_mean");
  double normal_sdev = Test_::test().getDouble("normal_sdev");
  util::NormalDistribution<double> z(N, normal_mean, normal_sdev, seed);
  std::vector<double> z_check = Test_::test().getDoubleVector("normal_double_answer");
  oops::Log::info() << "\nTesting oops::util::Random.h Gaussian Distribution: \n"
                    << z << std::endl;
  for (std::size_t jj = 0; jj < N; ++jj)
    EXPECT(oops::is_close(z[jj], z_check[jj],
                                                 tol * std::abs(z_check[jj])));

  /*! Test the shuffling algorithm */
  std::vector<int> shuffled_vector(N);
  std::iota(shuffled_vector.begin(), shuffled_vector.end(), 0);
  util::shuffle(shuffled_vector.begin(), shuffled_vector.end(), seed);
  std::vector<int> shuffled_vector_check = Test_::test().getIntVector("shuffle_answer");
  oops::Log::info() << "\nTesting oops::util::Random.h shuffle: \n"
                    << shuffled_vector << std::endl;
  EXPECT_EQUAL(shuffled_vector, shuffled_vector_check);
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

  EXPECT(test_uniform_real_f(&config) == 0);
  EXPECT(test_uniform_double_f(&config) == 0);
  EXPECT(test_uniform_int_f(&config) == 0);
  EXPECT(test_uniform_long_f(&config) == 0);
  EXPECT(test_normal_real_f(&config) == 0);
  EXPECT(test_normal_double_f(&config) == 0);
}

// -----------------------------------------------------------------------------

class Random : public oops::Test {
 public:
  using oops::Test::Test;
  virtual ~Random() {}

 private:
  std::string testid() const override {return "test::Random";}

  void register_tests() const override {
    std::vector<eckit::testing::Test>& ts = eckit::testing::specification();

    ts.emplace_back(CASE("util/Random/testCppRandom")
      { testCppRandom(); });
    ts.emplace_back(CASE("util/Random/testFortranRandom")
      { testFortranRandom(); });
  }

  void clear() const override {}
};

// -----------------------------------------------------------------------------
}  // namespace test

#endif  // TEST_UTIL_RANDOM_H_
