/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_FLOATCOMPARE_H_
#define OOPS_UTIL_FLOATCOMPARE_H_

#include <cmath>
#include <iomanip>
#include <limits>
#include <vector>

#include "eckit/types/FloatCompare.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

enum class TestVerbosity {
  SILENT,
  LOG_FAILURE_ONLY,
  LOG_SUCCESS_AND_FAILURE
};

/// \brief Tests two floating-point numbers for approximate equality using a relative tolerance,
/// optionally logging the test result.
///
/// \note The other overload (without the \p log_prefix_generator parameter) should be sufficient
/// for most applications.
///
/// \param a, b
///   Numbers to compare.
/// \param max_relative_difference
///   Maximum acceptable relative difference between \p a and \p b.
/// \param log_prefix_generator
///   A functor taking a single parameter of type std::ostream&. If the test result needs to be
///   logged, the functor will be called with a reference to the log stream passed to the argument.
///   This can be used to output a custom prefix (such as extra information about the compared
///   values) to the log stream.
/// \param verbosity
///   Determines whether the test result will be logged both on success and failure, only on
///   failure, or not at all.
/// \param max_ulps_diff
///   Maximum spacing between \p a and \p b in ULPs (where 1 ULP, or the unit of least precision,
///   is the spacing between two consecutive floating-point numbers) for which the function
///   returns true even if the relative difference between \p a and \p b is larger than
///   maximum_relative_difference. By default, 0.
///
/// \returns False if |a - b| > max(|a|, |b|) * max_relative_difference or if either \p a or \p b
/// is NaN or infinite; true otherwise. (wsmigaj: Not sure why false is returned for two infinities
/// of the same sign.)
template <typename T, typename LogPrefixGenerator>
bool is_close_relative(
    T a, T b, T max_relative_difference, int max_ulps_diff,
    const LogPrefixGenerator & log_prefix_generator,
    TestVerbosity verbosity = TestVerbosity::LOG_FAILURE_ONLY) {
  // if nan or inf values, always return false
  if (std::isnan(a) || std::isnan(b) || std::isinf(a) || std::isinf(b)) return false;

  // otherwise, create an absolute tolerance that is of the same type as a and b
  // (which is what is_approximately_equal wants) and call is_approximately_equal.
  T AbsA = fabs(a);
  T AbsB = fabs(b);
  T MaxAbs = (AbsA < AbsB ? AbsB : AbsA);
  // greater of AbsA, AbsB times max_relative_difference
  T EpsAB = MaxAbs * max_relative_difference;
  bool passed = eckit::types::is_approximately_equal(a, b, EpsAB, max_ulps_diff);
  std::size_t num_digits = std::numeric_limits<T>::max_digits10;
  if (passed) {
    if (verbosity == TestVerbosity::LOG_SUCCESS_AND_FAILURE) {
      log_prefix_generator(Log::info());
      Log::info() << "difference between " << std::setprecision(num_digits)
                  << a << " and " << b << " is " << fabs(a - b)
                  << ", i.e. no more than " << EpsAB << " (PASS)" << std::endl;
    }
  } else {
    if (verbosity != TestVerbosity::SILENT) {
      log_prefix_generator(Log::info());
      Log::info() << "difference between " << std::setprecision(num_digits)
                  << a << " and " << b << " is " << fabs(a - b)
                  << ", i.e. exceeds " << EpsAB << " (FAIL)" << std::endl;
    }
  }

  return passed;
}

/// \brief Tests two floating-point numbers for approximate equality using a relative tolerance,
/// optionally logging the test result.
///
/// \param a, b
///   Numbers to compare.
/// \param max_relative_difference
///   Maximum acceptable relative difference between \p a and \p b.
/// \param verbosity
///   Determines whether the test result will be logged both on success and failure, only on
///   failure, or not at all.
/// \param max_ulps_diff
///   Maximum spacing between \p a and \p b in ULPs (where 1 ULP, or the unit of least precision,
///   is the spacing between two consecutive floating-point numbers) for which the function
///   returns true even if the relative difference between \p a and \p b is larger than
///   maximum_relative_difference. By default, 0.
///
/// \returns False if |a - b| > max(|a|, |b|) * max_relative_difference or if either \p a or \p b
/// is NaN or infinite; true otherwise. (wsmigaj: Not sure why false is returned for two infinities
/// of the same sign.)
template <typename T>
bool is_close_relative(
    T a, T b, T max_relative_difference, int max_ulps_diff = 0,
    TestVerbosity verbosity = TestVerbosity::LOG_FAILURE_ONLY) {
  return is_close_relative(a,
                           b,
                           max_relative_difference,
                           max_ulps_diff,
                           [](std::ostream &) {},
                           verbosity);
}

/// \brief The same as is_close_relative. In new code, prefer the longer name as it's more explicit.
template <typename T>
bool is_close(T a, T b, T max_relative_difference, int max_ulps_diff = 0,
              TestVerbosity verbosity = TestVerbosity::LOG_FAILURE_ONLY) {
  return is_close_relative(a, b, max_relative_difference, max_ulps_diff, verbosity);
}

/// \brief Tests two vectors of floating-point numbers for approximate equality using a relative
/// tolerance, optionally logging the test result.
///
/// \param a, b
///   Vectors to compare.
/// \param max_relative_difference
///   Maximum acceptable relative difference between corresponding elements of \p a and \p b.
/// \param verbosity
///   Determines whether the test result will be logged both on success and failure, only on
///   failure, or not at all.
/// \param max_ulps_diff
///   Maximum spacing between \p a and \p b in ULPs (where 1 ULP, or the unit of least precision,
///   is the spacing between two consecutive vectors) for which the function
///   returns true even if the relative difference between \p a and \p b is larger than
///   maximum_relative_difference. By default, 0.
///
/// \returns False if \p and \p b have different lengths or if for any index i
/// |a[i] - b[i]| > max(|a[i]|, |b[i]|) * max_relative_difference or either \p a[i] or \p b[i]
/// is NaN or infinite; true otherwise. (wsmigaj: Not sure why false is returned for two infinities
/// of the same sign.)
template <typename T>
bool are_all_close_relative(
    const std::vector<T> &a, const std::vector<T> &b, T max_relative_difference,
    int max_ulps_diff = 0, TestVerbosity verbosity = TestVerbosity::LOG_FAILURE_ONLY) {
  if (a.size() != b.size()) {
    if (verbosity != TestVerbosity::SILENT) {
      Log::info() << "vector lengths (" << a.size() << ", " << b.size() << ") don't match (FAIL)"
                  << std::endl;
    }
    return false;
  }

  bool passed = true;
  for (size_t i = 0; i < a.size(); ++i) {
    if (!is_close_relative(
          a[i], b[i], max_relative_difference, max_ulps_diff,
          [i] (std::ostream &os) { os << "vector element #" << i << ": "; },
          verbosity))
      passed = false;
  }
  return passed;
}

/// \brief Tests two floating-point numbers for approximate equality using an absolute tolerance,
/// optionally logging the test result.
///
/// \note The other overload (without the \p log_prefix_generator parameter) should be sufficient
/// for most applications.
///
/// \param a, b
///   Numbers to compare.
/// \param max_absolute_difference
///   Maximum acceptable absolute difference between \p a and \p b.
/// \param log_prefix_generator
///   A functor taking a single parameter of type std::ostream&. If the test result needs to be
///   logged, the functor will be called with a reference to the log stream passed to the argument.
///   This can be used to output a custom prefix (such as extra information about the compared
///   values) to the log stream.
/// \param verbosity
///   Determines whether the test result will be logged both on success and failure, only on
///   failure, or not at all.
/// \param max_ulps_diff
///   Maximum spacing between \p a and \p b in ULPs (where 1 ULP, or the unit of least precision,
///   is the spacing between two consecutive floating-point numbers) for which the function
///   returns true even if the relative difference between \p a and \p b is larger than
///   maximum_relative_difference. By default, 0.
///
/// \returns False if |a - b| > max_absolute_difference or if either \p a or \p b is NaN or
/// infinite; true otherwise. (wsmigaj: Not sure why false is returned for two infinities of the
/// same sign.)
template <typename T, typename LogPrefixGenerator>
bool is_close_absolute(
    T a, T b, T max_absolute_difference, int max_ulps_diff,
    const LogPrefixGenerator & log_prefix_generator,
    TestVerbosity verbosity = TestVerbosity::LOG_FAILURE_ONLY) {
  bool passed = eckit::types::is_approximately_equal(a, b, max_absolute_difference, max_ulps_diff);
  std::size_t num_digits = std::numeric_limits<T>::max_digits10;
  if (passed) {
    if (verbosity == TestVerbosity::LOG_SUCCESS_AND_FAILURE) {
      log_prefix_generator(Log::info());
      Log::info() << "difference between " << std::setprecision(num_digits)
                  << a << " and " << b << " is " << fabs(a - b)
                  << ", i.e. no more than " << max_absolute_difference << " (PASS)" << std::endl;
    }
  } else {
    if (verbosity != TestVerbosity::SILENT) {
      log_prefix_generator(Log::info());
      Log::info() << "difference between " << std::setprecision(num_digits)
                  << a << " and " << b << " is " << fabs(a - b)
                  << ", i.e. exceeds " << max_absolute_difference << " (FAIL)" << std::endl;
    }
  }

  return passed;
}

/// \brief Tests two floating-point numbers for approximate equality using an absolute tolerance,
/// optionally logging the test result.
///
/// \param a, b
///   Numbers to compare.
/// \param max_absolute_difference
///   Maximum acceptable absolute difference between \p a and \p b.
/// \param verbosity
///   Determines whether the test result will be logged both on success and failure, only on
///   failure, or not at all.
/// \param max_ulps_diff
///   Maximum spacing between \p a and \p b in ULPs (where 1 ULP, or the unit of least precision,
///   is the spacing between two consecutive floating-point numbers) for which the function
///   returns true even if the relative difference between \p a and \p b is larger than
///   maximum_relative_difference. By default, 0.
///
/// \returns False if |a - b| > max_absolute_difference or if either \p a or \p b is NaN or
/// infinite; true otherwise. (wsmigaj: Not sure why false is returned for two infinities of the
/// same sign.)
template <typename T>
bool is_close_absolute(
    T a, T b, T max_relative_difference, int max_ulps_diff = 0,
    TestVerbosity verbosity = TestVerbosity::LOG_FAILURE_ONLY) {
  return is_close_absolute(a,
                           b,
                           max_relative_difference,
                           max_ulps_diff,
                           [](std::ostream &) {},
                           verbosity);
}

/// \brief Tests two vectors of floating-point numbers for approximate equality using an absolute
/// tolerance, optionally logging the test result.
///
/// \param a, b
///   Vectors to compare.
/// \param max_relative_difference
///   Maximum acceptable absolute difference between corresponding elements of \p a and \p b.
/// \param verbosity
///   Determines whether the test result will be logged both on success and failure, only on
///   failure, or not at all.
/// \param max_ulps_diff
///   Maximum spacing between \p a and \p b in ULPs (where 1 ULP, or the unit of least precision,
///   is the spacing between two consecutive vectors) for which the function
///   returns true even if the relative difference between \p a and \p b is larger than
///   maximum_relative_difference. By default, 0.
///
/// \returns False if \p and \p b have different lengths or if for any index i |a[i] - b[i]| >
/// max_absolute_difference or either \p a[i] or \p b[i] is NaN or infinite; true otherwise.
/// (wsmigaj: Not sure why false is returned for two infinities of the same sign.)
template <typename T>
bool are_all_close_absolute(
    const std::vector<T> &a, const std::vector<T> &b, T max_absolute_difference,
    int max_ulps_diff = 0, TestVerbosity verbosity = TestVerbosity::LOG_FAILURE_ONLY) {
  if (a.size() != b.size()) {
    if (verbosity != TestVerbosity::SILENT) {
      Log::info() << "vector lengths (" << a.size() << ", " << b.size() << ") don't match (FAIL)"
                  << std::endl;
    }
    return false;
  }

  bool passed = true;
  for (size_t i = 0; i < a.size(); ++i) {
    if (!is_close_absolute(
          a[i], b[i], max_absolute_difference, max_ulps_diff,
          [i] (std::ostream &os) { os << "vector element #" << i << ": "; },
          verbosity))
      passed = false;
  }
  return passed;
}

}  // namespace oops
#endif  // OOPS_UTIL_FLOATCOMPARE_H_
