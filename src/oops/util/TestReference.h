/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_TESTREFERENCE_H_
#define OOPS_UTIL_TESTREFERENCE_H_

#include <exception>
#include <fstream>
#include <iostream>
#include <string>

#include "eckit/config/LocalConfiguration.h"


namespace oops {

// -----------------------------------------------------------------------------
class TestReference {
 public:
  TestReference() = default;  // default constructor

  void testReferenceFinalise(std::stringstream &);
  void Initialise(eckit::LocalConfiguration &);

 private:
  std::ifstream readRefFile();
  void testCompare(std::stringstream &, std::ifstream &, float, int);
  bool initCheck_ = false;

  std::string refFileName_;
  std::string outputFileYaml_;
  std::string refFileYaml_;
  std::string testFileYaml_;
  double tolFloat_;
  int tolInt_;
};
// -----------------------------------------------------------------------------
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

// Textual mismatch error condition
class TestReferenceTextMismatchError : public TestReferenceError
{
 public:
  TestReferenceTextMismatchError(int line_num,
  const std::string &test_line, const std::string &ref_line);
};

// Integer mismatch error condition
class TestReferenceIntegerMismatchError : public TestReferenceError {
 public:
  using NumT = int64_t;
    TestReferenceIntegerMismatchError(int line_num,
          NumT test_val, NumT ref_val, NumT tolerance,
          const std::string &test_line, const std::string &ref_line);
};

// Floating-point mismatch error condition
class TestReferenceFloatMismatchError : public TestReferenceError {
 public:
    using NumT = double;
    TestReferenceFloatMismatchError(int line_num,
          NumT test_val, NumT ref_val, NumT tolerance,
          const std::string &test_line, const std::string &ref_line);
};


}  // namespace oops

#endif  // OOPS_UTIL_TESTREFERENCE_H_
