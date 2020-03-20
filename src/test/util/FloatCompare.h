/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_FLOATCOMPARE_H_
#define TEST_UTIL_FLOATCOMPARE_H_

#include <limits>
#include <map>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/parallel/mpi/mpi.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"

namespace test {

template <typename T>
void testIsRelativeDifferenceAtMost(oops::TestVerbosity verbosity)
{
  const T nan = std::numeric_limits<T>::quiet_NaN();
  const T inf = std::numeric_limits<T>::infinity();

  // Positive numbers
  EXPECT(oops::is_close_relative(T(2.0), T(4.0), T(0.5), verbosity));
  EXPECT(oops::is_close_relative(T(4.0), T(2.0), T(0.5), verbosity));
  EXPECT_NOT(oops::is_close_relative(T(2.0), T(4.0), T(0.49), verbosity));
  EXPECT_NOT(oops::is_close_relative(T(4.0), T(2.0), T(0.49), verbosity));

  // Negative numbers
  EXPECT(oops::is_close_relative(T(-2.0), T(-4.0), T(0.5), verbosity));
  EXPECT(oops::is_close_relative(T(-4.0), T(-2.0), T(0.5), verbosity));
  EXPECT_NOT(oops::is_close_relative(T(-2.0), T(-4.0), T(0.49), verbosity));
  EXPECT_NOT(oops::is_close_relative(T(-4.0), T(-2.0), T(0.49), verbosity));

  // NaNs
  EXPECT_NOT(oops::is_close_relative(nan, T(1.0), T(0.1), verbosity));
  EXPECT_NOT(oops::is_close_relative(T(1.0), nan, T(0.1), verbosity));
  EXPECT_NOT(oops::is_close_relative(nan, nan, T(0.1), verbosity));

  // Infinities
  EXPECT_NOT(oops::is_close_relative(inf, T(1.0), T(0.1), verbosity));
  EXPECT_NOT(oops::is_close_relative(T(1.0), inf, T(0.1), verbosity));
}

template <typename T>
void testIsRelativeDifferenceAtMost()
{
  oops::Log::info() << "In the following, neither successes nor failures should be logged\n";
  testIsRelativeDifferenceAtMost<T>(oops::TestVerbosity::SILENT);
  oops::Log::info() << "In the following, only failures should be logged\n";
  testIsRelativeDifferenceAtMost<T>(oops::TestVerbosity::LOG_FAILURE_ONLY);
  oops::Log::info() << "In the following, both successes and failures should be logged\n";
  testIsRelativeDifferenceAtMost<T>(oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE);
}

CASE("util/FloatComparisons/is_close_relative/float") {
  testIsRelativeDifferenceAtMost<float>();
}

CASE("util/FloatComparisons/is_close_relative/double") {
  testIsRelativeDifferenceAtMost<double>();
}

template <typename T>
void testIsAbsoluteDifferenceAtMost(oops::TestVerbosity verbosity)
{
  const T nan = std::numeric_limits<T>::quiet_NaN();
  const T inf = std::numeric_limits<T>::infinity();

  // Positive numbers
  EXPECT(oops::is_close_absolute(T(2.0), T(4.0), T(2.0), verbosity));
  EXPECT(oops::is_close_absolute(T(4.0), T(2.0), T(2.0), verbosity));
  EXPECT_NOT(oops::is_close_absolute(T(2.0), T(4.0), T(1.99), verbosity));
  EXPECT_NOT(oops::is_close_absolute(T(4.0), T(2.0), T(1.99), verbosity));

  // Negative numbers
  EXPECT(oops::is_close_absolute(T(-2.0), T(-4.0), T(2.0), verbosity));
  EXPECT(oops::is_close_absolute(T(-4.0), T(-2.0), T(2.0), verbosity));
  EXPECT_NOT(oops::is_close_absolute(T(-2.0), T(-4.0), T(1.99), verbosity));
  EXPECT_NOT(oops::is_close_absolute(T(-4.0), T(-2.0), T(1.99), verbosity));

  // NaNs
  EXPECT_NOT(oops::is_close_absolute(nan, T(1.0), T(0.1), verbosity));
  EXPECT_NOT(oops::is_close_absolute(T(1.0), nan, T(0.1), verbosity));
  EXPECT_NOT(oops::is_close_absolute(nan, nan, T(0.1), verbosity));

  // Infinities
  EXPECT_NOT(oops::is_close_absolute(inf, T(1.0), T(0.1), verbosity));
  EXPECT_NOT(oops::is_close_absolute(T(1.0), inf, T(0.1), verbosity));
}

template <typename T>
void testIsAbsoluteDifferenceAtMost()
{
  oops::Log::info() << "In the following, neither successes nor failures should be logged\n";
  testIsAbsoluteDifferenceAtMost<T>(oops::TestVerbosity::SILENT);
  oops::Log::info() << "In the following, only failures should be logged\n";
  testIsAbsoluteDifferenceAtMost<T>(oops::TestVerbosity::LOG_FAILURE_ONLY);
  oops::Log::info() << "In the following, both successes and failures should be logged\n";
  testIsAbsoluteDifferenceAtMost<T>(oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE);
}

CASE("util/FloatComparisons/is_close_absolute/float") {
  testIsAbsoluteDifferenceAtMost<float>();
}

CASE("util/FloatComparisons/is_close_absolute/double") {
  testIsAbsoluteDifferenceAtMost<double>();
}

template <typename T>
void testAreAllRelativeDifferencesAtMost(oops::TestVerbosity verbosity)
{
  // Same lengths
  EXPECT(oops::are_all_close_relative(
           std::vector<T>{}, std::vector<T>{}, T(0.5), verbosity));
  EXPECT(oops::are_all_close_relative(
           std::vector<T>{T(2.0)}, std::vector<T>{T(4.0)}, T(0.5), verbosity));
  EXPECT(oops::are_all_close_relative(
           std::vector<T>{T(2.0), T(-2.0)}, std::vector<T>{T(4.0), T(-4.0)}, T(0.5), verbosity));
  EXPECT_NOT(oops::are_all_close_relative(
               std::vector<T>{T(1.0), T(-2.0)}, std::vector<T>{T(4.0), T(-4.0)},
               T(0.5), verbosity));
  EXPECT_NOT(oops::are_all_close_relative(
               std::vector<T>{T(2.0), T(-1.0)}, std::vector<T>{T(4.0), T(-4.0)},
               T(0.5), verbosity));

  // Different lengths
  EXPECT_NOT(oops::are_all_close_relative(
           std::vector<T>{}, std::vector<T>{1.0}, T(0.5), verbosity));
  EXPECT_NOT(oops::are_all_close_relative(
           std::vector<T>{1.0}, std::vector<T>{}, T(0.5), verbosity));
}

template <typename T>
void testAreAllRelativeDifferencesAtMost()
{
  oops::Log::info() << "In the following, neither successes nor failures should be logged\n";
  testAreAllRelativeDifferencesAtMost<T>(oops::TestVerbosity::SILENT);
  oops::Log::info() << "In the following, only failures should be logged\n";
  testAreAllRelativeDifferencesAtMost<T>(oops::TestVerbosity::LOG_FAILURE_ONLY);
  oops::Log::info() << "In the following, both successes and failures should be logged\n";
  testAreAllRelativeDifferencesAtMost<T>(oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE);
}

CASE("util/FloatComparisons/are_all_close_relative/float") {
  testAreAllRelativeDifferencesAtMost<float>();
}

CASE("util/FloatComparisons/are_all_close_relative/double") {
  testAreAllRelativeDifferencesAtMost<double>();
}

template <typename T>
void testAreAllAbsoluteDifferencesAtMost(oops::TestVerbosity verbosity)
{
  // Same lengths
  EXPECT(oops::are_all_close_absolute(
           std::vector<T>{}, std::vector<T>{}, T(2.0), verbosity));
  EXPECT(oops::are_all_close_absolute(
           std::vector<T>{T(2.0)}, std::vector<T>{T(4.0)}, T(2.0), verbosity));
  EXPECT(oops::are_all_close_absolute(
           std::vector<T>{T(2.0), T(-2.0)}, std::vector<T>{T(4.0), T(-4.0)}, T(2.0), verbosity));
  EXPECT_NOT(oops::are_all_close_absolute(
               std::vector<T>{T(1.0), T(-2.0)}, std::vector<T>{T(4.0), T(-4.0)},
               T(2.0), verbosity));
  EXPECT_NOT(oops::are_all_close_absolute(
               std::vector<T>{T(2.0), T(-1.0)}, std::vector<T>{T(4.0), T(-4.0)},
               T(2.0), verbosity));

  // Different lengths
  EXPECT_NOT(oops::are_all_close_absolute(
           std::vector<T>{}, std::vector<T>{1.0}, T(2.0), verbosity));
  EXPECT_NOT(oops::are_all_close_absolute(
           std::vector<T>{1.0}, std::vector<T>{}, T(2.0), verbosity));
}

template <typename T>
void testAreAllAbsoluteDifferencesAtMost()
{
  oops::Log::info() << "In the following, neither successes nor failures should be logged\n";
  testAreAllAbsoluteDifferencesAtMost<T>(oops::TestVerbosity::SILENT);
  oops::Log::info() << "In the following, only failures should be logged\n";
  testAreAllAbsoluteDifferencesAtMost<T>(oops::TestVerbosity::LOG_FAILURE_ONLY);
  oops::Log::info() << "In the following, both successes and failures should be logged\n";
  testAreAllAbsoluteDifferencesAtMost<T>(oops::TestVerbosity::LOG_SUCCESS_AND_FAILURE);
}

CASE("util/FloatComparisons/are_all_close_absolute/float") {
  testAreAllAbsoluteDifferencesAtMost<float>();
}

CASE("util/FloatComparisons/are_all_close_absolute/double") {
  testAreAllAbsoluteDifferencesAtMost<double>();
}

class FloatCompare : public oops::Test {
 private:
  std::string testid() const override {return "test::FloatCompare";}

  void register_tests() const override {}
};

}  // namespace test

#endif  // TEST_UTIL_FLOATCOMPARE_H_
