/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/util/TestReference.h"

#include <algorithm>
#include <cctype>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>
#include <typeinfo>
#include <utility>
#include <vector>

#if USE_BOOST_REGEX
#include <boost/regex.hpp>
#define REGEX_NAMESPACE boost
#else
#include <regex>
#define REGEX_NAMESPACE std
#endif

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/mpi/mpi.h"
#include "oops/util/LibOOPS.h"
#include "oops/util/Logger.h"

namespace oops {

/** Local helper routines to parse numbers without using std::regex facilities
 */
namespace {

// -----------------------------------------------------------------------------
  /** Determine if string is a valid integer representation
   */
  bool is_integer_repr(const std::string &s)
  {
    auto begin = s.begin();
    // Drop initial +/- char
    if (*begin == '-' || *begin == '+') {
      begin++;
      if (begin == s.end()) return false;
    }
    // All remaining chars should be digits
    return std::all_of(begin, s.end(), [](unsigned char c){return std::isdigit(c);} );
  }

// -----------------------------------------------------------------------------
  /** Return a vector of all non-overlapping substrings that can be parsed as valid numbers
   */
  std::vector<std::string> parse_numbers(const std::string &s)
  {
    std::vector<std::string> numbers;
    std::string init_char = "+-0123456789";  // A number must start with these chars
    auto start = s.find_first_of(init_char);
    while (start < s.size()) {
      // Number must have a digit as first or second char
      if (!std::isdigit(s[start])) {
        if (start == s.size()-1) break;
        if (!std::isdigit(s[start+1])) {
          start = s.find_first_of(init_char, start+1);
          continue;
        }
      }
      // Sequence of chars will parse as a number, pos is the extent successfully converted
      size_t pos = 0;
      std::stod(s.substr(start), &pos);
      numbers.push_back(s.substr(start, pos));
      start = s.find_first_of(init_char, start+pos);
    }
    return numbers;
  }

// -----------------------------------------------------------------------------
  /** Trim a string (remove leading and trailing whitespaces)
   */
  std::string trim(const std::string& line)
  {
    const char* whitespace = " \t\v\r\n";
    const auto start = line.find_first_not_of(whitespace);
    const auto end = line.find_last_not_of(whitespace);
    return start == end ? std::string() : line.substr(start, end - start + 1);
  }

}  // anonymous namespace

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

void TestReference::initialise(const eckit::LocalConfiguration &conf)
{
  initCheck_ = true;
  std::string mpi_pattern, mpi_size;
  if (conf.has("mpi pattern")) {
    mpi_pattern = conf.getString("mpi pattern");
    mpi_size = std::to_string(oops::mpi::world().size());
  }

  refFile_ = conf.getString("reference filename");
  if (conf.has("mpi pattern")) {
    refFile_  = REGEX_NAMESPACE::regex_replace(refFile_,
                                               REGEX_NAMESPACE::regex(mpi_pattern),
                                               mpi_size);
  }
  oops::Log::info() << "[TestReference] Comparing to reference file: " << refFile_ << std::endl;

  tolRelativeFloat_ = conf.getDouble("float relative tolerance", 1.0e-6);
  tolAbsoluteFloat_ = conf.getDouble("float absolute tolerance", 0.0);
  oops::Log::info() << "Relative float tolerance for tests : " << tolRelativeFloat_
                     << std::endl;
  oops::Log::info() << "Absolute float tolerance for tests : " << tolAbsoluteFloat_
                     << std::endl;

  tolInt_ = conf.getInt("integer tolerance", 0);

  if (conf.has("log output filename")) {
    outputFile_ = conf.getString("log output filename");
    if (conf.has("mpi pattern")) {
      outputFile_  = REGEX_NAMESPACE::regex_replace(outputFile_,
                                                    REGEX_NAMESPACE::regex(mpi_pattern),
                                                    mpi_size);
    }
    LibOOPS::instance().teeOutput(outputFile_);
    oops::Log::info() << "[TestReference] Saving Log output to: " << outputFile_ << std::endl;
  }

  if (conf.has("test output filename")) {
    testFile_ = conf.getString("test output filename");
    if (conf.has("mpi pattern")) {
      testFile_  = REGEX_NAMESPACE::regex_replace(testFile_,
                                                  REGEX_NAMESPACE::regex(mpi_pattern),
                                                  mpi_size);
    }
    oops::Log::info() << "[TestReference] Saving Test output to: " << testFile_ << std::endl;
  }
}

// -----------------------------------------------------------------------------

void TestReference::compare(const std::string & test,
                            const std::string & ref,
                            FloatT tolRelativeFloat,
                            FloatT tolAbsoluteFloat,
                            IntT tolInt)
{
  int lineCounter = 0;
  std::string testPrefix = "Test     : ";  // Test prefix string to remove
  std::stringstream testStream(test);
  std::stringstream refStream(ref);
  std::string refLine, testLine;

  while (std::getline(refStream, refLine)) {
    lineCounter++;

    // If reference line is empty, go to the next line
    if (refLine.empty()) continue;

    testLine = "";
    // Find non-empty line in the test output
    while (testLine.empty()) {
      std::getline(testStream, testLine);
      // Check that a corresponding test line exists
      if (testStream.fail()) {
        throw TestReferenceMissingTestLineError(lineCounter, refLine);
      }
    }

    // Remove test prefix from refLine if present
    if (testPrefix == refLine.substr(0, testPrefix.size())) {
      refLine.erase(0, testPrefix.size());
    }

    if (trim(refLine) == trim(testLine)) continue;

    // If strings don't match - parse all numbers in the test and reference strings
    auto refNums = parse_numbers(refLine);
    auto testNums = parse_numbers(testLine);

    if (refNums.empty() || testNums.empty() || refNums.size() != testNums.size()) {
      // No numbers parsed or differing numbers parsed in each line.
      throw TestReferenceTextMismatchError(lineCounter, testLine, refLine);
    }

    // Check if each pair of numbers are within tolerance
    for (size_t i = 0; i < refNums.size(); i++) {
      if (is_integer_repr(refNums[i]) && is_integer_repr(testNums[i])) {
        // both test and reference appear to be integers
        auto ref_int = stol(refNums[i]);
        auto test_int = stol(testNums[i]);
        if (std::abs(ref_int - test_int) > tolInt)
          throw TestReferenceIntegerMismatchError(lineCounter, test_int,
            ref_int, std::abs(ref_int - test_int), tolInt, testLine, refLine);
      } else {
        FloatT ref_float = stod(refNums[i]);
        FloatT test_float = stod(testNums[i]);
        FloatT diff = std::abs(ref_float - test_float);
        if (diff > tolRelativeFloat * std::abs(0.5 * (ref_float + test_float)) &&
            diff > tolAbsoluteFloat) {
          throw TestReferenceFloatMismatchError(lineCounter, test_float, ref_float, diff,
            tolRelativeFloat * std::abs(0.5 * (ref_float + test_float)),
            tolAbsoluteFloat, testLine, refLine);
        }
      }
    }
  }

  // Check there are no more remaining test lines to process
  testLine = "";
  // Read until either reach end of file or find non-empty line in the test output
  while (testLine.empty() && (!testStream.fail())) {
    std::getline(testStream, testLine);
  }
  if (!testLine.empty()) {
    throw TestReferenceMissingReferenceLineError(lineCounter, testLine);
  }
}

// -----------------------------------------------------------------------------

void TestReference::finalise(const std::string & testStr)
{
  if (!initCheck_) return;

  // Read reference file to string
  std::ifstream refFileIn(refFile_);
  if (refFileIn.fail()) throw eckit::CantOpenFile(refFile_);
  std::string refStr(std::istreambuf_iterator<char>{refFileIn}, std::istreambuf_iterator<char>{});
  refFileIn.close();

  if (!testFile_.empty()) {
    // Write test string to file
    std::ofstream testFileOut(testFile_);
    if (!testFileOut) throw eckit::CantOpenFile(testFile_);
    testFileOut << testStr;
    testFileOut.close();
  }

  compare(testStr, refStr, tolRelativeFloat_, tolAbsoluteFloat_, tolInt_);
  oops::Log::info() << "[TestReference] Comparison is done" << std::endl;
}

// -----------------------------------------------------------------------------
// -----------------------------------------------------------------------------

TestReferenceMissingReferenceLineError::TestReferenceMissingReferenceLineError(int line_num,
  const std::string &test_line)
{
  std::ostringstream os;
  os << "TestReference: Missing reference file line corresponding to test output Line#:" << line_num
     << "\nTest line: '" << test_line << "'";
  what_ = os.str();
}

// -----------------------------------------------------------------------------

TestReferenceMissingTestLineError::TestReferenceMissingTestLineError(int line_num,
  const std::string &ref_line)
{
  std::ostringstream os;
  os << "TestReference: Missing test output line corresponding to reference file Line#:" << line_num
     << "\nRef line: '" << ref_line << "'";
  what_ = os.str();
}

// -----------------------------------------------------------------------------

TestReferenceTextMismatchError::TestReferenceTextMismatchError(int line_num,
  const std::string &test_line, const std::string &ref_line)
{
  std::ostringstream os;
  os << "Test reference Text mismatch @ Line:" << line_num
     << "\nTest: '" << test_line << "'\nRef:  '" << ref_line << "'";
  what_ = os.str();
}

// -----------------------------------------------------------------------------

TestReferenceIntegerMismatchError::TestReferenceIntegerMismatchError(int line_num,
                                  NumT test_val, NumT ref_val, NumT diff, NumT tolerance,
                                  const std::string &test_line,
                                  const std::string &ref_line)
{
  std::ostringstream os;
  os << "Test reference Integer mismatch @ Line:" << line_num << "\n"
     << "Test Val : " << test_val << "\n"
     << "Ref  Val : " << ref_val << "\n"
     << "Delta    : " << diff << "\n"
     << "Tolerance: " << tolerance << "\n"
     << "Test Line: '" << test_line << "'\n"
     << "Ref Line : '" << ref_line << "'";
  what_ = os.str();
}

// -----------------------------------------------------------------------------

TestReferenceFloatMismatchError::TestReferenceFloatMismatchError(int line_num,
                                  NumT test_val, NumT ref_val, NumT diff,
                                  NumT relativeTolerance,
                                  NumT absoluteTolerance,
                                  const std::string &test_line,
                                  const std::string &ref_line)
{
  std::ostringstream os;
  os << "Test reference Float mismatch @ Line:" << line_num <<"\n"
     << std::setprecision(std::numeric_limits<NumT>::digits10 + 1)
     << std::scientific
     << "Test Val : " << test_val << "\n"
     << "Ref  Val : " << ref_val << "\n"
     << "Delta    : "
     << diff
     << "\n"
     << "Relative tolerance: " << relativeTolerance << "\n"
     << "Absolute tolerance: " << absoluteTolerance << "\n"
     << "Test Line: '" << test_line << "'\n"
     << "Ref Line : '" << ref_line << "'";
  what_ = os.str();
}

}  // namespace oops
