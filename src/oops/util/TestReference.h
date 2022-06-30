/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_TESTREFERENCE_H_
#define OOPS_UTIL_TESTREFERENCE_H_

#include <cstdint>
#include <exception>
#include <string>

namespace eckit {
  class LocalConfiguration;
}

namespace oops {

class TestReference {
 public:
  using FloatT = double;
  using IntT = int64_t;

  void finalise(const std::string &testStr);
  void initialise(const eckit::LocalConfiguration &conf);

  static void compare(const std::string &, const std::string &,
                      FloatT, FloatT, IntT);

 private:
  bool initCheck_ = false;

  std::string refFile_;
  std::string outputFile_;
  std::string testFile_;
  FloatT tolRelativeFloat_;
  FloatT tolAbsoluteFloat_;
  IntT tolInt_;
};

// Base class for all TestReference errors
class TestReferenceError : public std::exception {
 public:
  const char* what() const noexcept override
  {
    return what_.c_str();
  }

 protected:
  std::string what_;  // error message to print
};

// Error: One or more reference lines are missing
class TestReferenceMissingReferenceLineError : public TestReferenceError
{
 public:
  TestReferenceMissingReferenceLineError(int line_num,
  const std::string &test_line);
};

// Error: One or more test lines are missing
class TestReferenceMissingTestLineError : public TestReferenceError
{
 public:
  TestReferenceMissingTestLineError(int line_num,
  const std::string &ref_line);
};

// Error: Test and reference lines don't match
class TestReferenceTextMismatchError : public TestReferenceError
{
 public:
  TestReferenceTextMismatchError(int line_num,
  const std::string &test_line, const std::string &ref_line);
};

// Error: Test and reference lines parsed integer representations not within tolerance
class TestReferenceIntegerMismatchError : public TestReferenceError {
 public:
  using NumT = TestReference::IntT;
    TestReferenceIntegerMismatchError(int line_num,
          NumT test_val, NumT ref_val, NumT diff,  NumT tolerance,
          const std::string &test_line, const std::string &ref_line);
};

// Error: Test and reference lines parsed floating-point representations not within tolerance
class TestReferenceFloatMismatchError : public TestReferenceError {
 public:
    using NumT = TestReference::FloatT;
    TestReferenceFloatMismatchError(int line_num,
          NumT test_val, NumT ref_rel_val, NumT ref_abs_val, NumT diff, NumT tolerance,
          const std::string &test_line, const std::string &ref_line);
};


}  // namespace oops

#endif  // OOPS_UTIL_TESTREFERENCE_H_
