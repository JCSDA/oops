/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#include <cmath>
#include <iomanip>
#include <iostream>
#include <limits>
#include <regex>
#include <sstream>
#include <typeinfo>
#include <vector>

#include "eckit/config/YAMLConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/util/LibOOPS.h"
#include "oops/util/Logger.h"
#include "oops/util/TestReference.h"


namespace oops {

// -----------------------------------------------------------------------------

void TestReference::Initialise(eckit::LocalConfiguration &testConf)
{
  initCheck_ = true;
  tolFloat_ = testConf.getFloat("float relative tolerance", 0.0);
  tolInt_ = testConf.getInt("integer tolerance", 0);
  refFileYaml_ = testConf.getString("reference filename");
  if (testConf.has("log output filename")) {
    outputFileYaml_ = testConf.getString("log output filename");
    LibOOPS::instance().teeOutput(outputFileYaml_);
  if (testConf.has("test output filename")) {
    testFileYaml_ = testConf.getString("test output filename");
    oops::Log::info() << "save test output to " << testFileYaml_ << std::endl;
  }
  }
}

// -----------------------------------------------------------------------------
// Reads reference file (refFileYaml_) for tests with comparing step
// refFileYaml_ is read from Yaml using setYaml method
std::ifstream TestReference::readRefFile()
{
  std::ifstream refFile(refFileYaml_);
  if (refFile.fail()) {
    throw eckit::CantOpenFile(refFileName_ + " does not exist");
  }
  return refFile;
}
// -----------------------------------------------------------------------------
void TestReference::testCompare(std::stringstream & test,
                           std::ifstream & ref,
                           float tolFloat, int tolInt)
{
  std::regex regNumber("[-+]?[0-9]*[.]?[0-9]+([eE][-+]?[0-9]+)?");
  std::regex regOnlyFloat("[-+]?[0-9]*[.][0-9]+([eE][-+]?[0-9]+)?");
  int lineCounter = 0;
  std::string refLine, testLine;
  while (std::getline(ref, refLine)) {
    std::getline(test, testLine);
    std::string testPrefix = "Test     : ";

    // if test prefix is in reference file add it to the testLine
    // to match the format
    // This can be removed once all reference files are updated and
    // test prefix is removed from them.
    if (!refLine.find(testPrefix)) {
      testLine = testPrefix+testLine;
    }

    lineCounter++;
    int compareRes = refLine.compare(testLine);

    // if strings don't match investigate more
    if (compareRes != 0) {
       // find all numbers in the two test and reference strings
       std::vector<std::string> numMatch_refLine(
         std::sregex_token_iterator(refLine.begin(), refLine.end(), regNumber),
         std::sregex_token_iterator());
       std::vector<std::string> numMatch_testLine(
         std::sregex_token_iterator(testLine.begin(), testLine.end(), regNumber),
         std::sregex_token_iterator());

       // How many numbers are in each string and do they match?
       if (numMatch_refLine.size() != numMatch_testLine.size()) {
         throw TestReferenceTextMismatchError(lineCounter, testLine, refLine);
       } else {
         for (auto i = 0; i < numMatch_refLine.size(); i++) {
           std::smatch matchCheck;
           //  separate float and int
           std::regex_search(numMatch_refLine[i], matchCheck, regOnlyFloat);

           if (matchCheck.empty()) {
             // number is int
             int ref_int = stoi(numMatch_refLine[i]);
             int test_int = stoi(numMatch_testLine[i]);
             int abs_dif = std::abs(ref_int - test_int);
             if (abs_dif > tolInt) {
               throw TestReferenceIntegerMismatchError(lineCounter, test_int,
                 ref_int, tolInt, testLine, refLine);
             }
           } else {
             // number is float
             double ref_float = stod(numMatch_refLine[i]);
             double test_float = stod(numMatch_testLine[i]);
             double rel_dir = std::abs((ref_float - test_float)/(0.5 * (ref_float + test_float)));
             if (rel_dir > tolFloat) {
               throw TestReferenceFloatMismatchError(lineCounter, test_float,
                 ref_float, tolFloat, testLine, refLine);
             }
           }
         }
       }
    }
  }
}
// -----------------------------------------------------------------------------
void TestReference::testReferenceFinalise(std::stringstream & testStream_)
{
  if (!initCheck_) {
    return;
  }
  auto refStream = readRefFile();
  if (!testFileYaml_.empty()) {
    std::ofstream testFileOut(testFileYaml_);
    testFileOut << testStream_.str();
    testFileOut.close();
  }
  testCompare(testStream_, refStream, tolFloat_, tolInt_);
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
                                  NumT test_val, NumT ref_val, NumT tolerance,
                                  const std::string &test_line,
                                  const std::string &ref_line)
{
    std::ostringstream os;
    os << "Test reference Integer mismatch @ Line:" << line_num << "\n"
       << "Test Val : " << test_val << "\n"
       << "Ref  Val : " << ref_val << "\n"
       << "Delta    : " << std::abs(test_val-ref_val) << "\n"
       << "Tolerance: " << tolerance << "\n"
       << "Test Line: '" << test_line << "'\n"
       << "Ref Line : '" << ref_line << "'";
    what_ = os.str();
}
// -----------------------------------------------------------------------------
TestReferenceFloatMismatchError::TestReferenceFloatMismatchError(int line_num,
                                  NumT test_val, NumT ref_val, NumT tolerance,
                                  const std::string &test_line, const std::string &ref_line)
{
    std::ostringstream os;
    os << "Test reference Float mismatch @ Line:" << line_num <<"\n"
       << std::setprecision(std::numeric_limits<NumT>::digits10 + 1)
       << std::scientific
       << "Test Val : " << test_val << "\n"
       << "Ref  Val : " << ref_val << "\n"
       << "Delta    : " << std::abs(test_val-ref_val) << "\n"
       << "Tolerance: " << tolerance << "\n"
       << "Test Line: '" << test_line << "'\n"
       << "Ref Line : '" << ref_line << "'";
    what_ = os.str();
}


}  // namespace oops
