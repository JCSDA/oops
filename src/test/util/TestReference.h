/*
 * (C) Copyright 2017-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_TESTREFERENCE_H_
#define TEST_UTIL_TESTREFERENCE_H_

#include <string>

#include "oops/runs/Test.h"
#include "oops/util/TestReference.h"

namespace test {

CASE("util/TestReference") {
  std::string test1 = "Line 1: +0123456789\n";
  std::string test2 = "Line 1: +0123456789\nLine 2: +99.98764321e-3\n";
  std::string test2_1 = "\n\nLine 1: +0123456789\n\n\nLine 2: +99.98764321e-3\n\n\n\n";
  std::string test3 = "-ABC-XYZ-\n";
  std::string test4 = "-99-+EE-++44E0\n";  // Parsed as {"-99", "+44E0"}
  std::string test5 = "0.123456789";

  std::string good_ref1 = "   Line 1: 123456789\t\n";
  std::string good_ref2 = "Line 1: 123456789\nLine 2: 9.9987643212e-2\n\n";
  std::string good_ref4 = "-9.90000E1-+EE-+44\n";
  std::string good_ref5 = "0.1234567885";

  std::string bad_ref1 = "Line 1: -123456789\n";
  std::string bad_ref2_1 = "Line 1: 123456789\nLine 2: 9.9987643e-2\n";
  std::string bad_ref2_2 = "Line 1: 123456789\nLine 2: 9.9987643212e-2 42\n";
  std::string bad_ref3 = "-ABC+XYZ-\n";
  std::string bad_ref4_1 = "-99-+EE-+-43\n";
  std::string bad_ref4_2 = "-90-+EE-+-43\n";
  std::string bad_ref4_3 = "-99-+EE-++44EE0\n";  // Parsed as {"-99", "+44", "0"}
  std::string bad_ref5 = "0.123456787";

  oops::TestReference::IntT iTol = 0;
  oops::TestReference::FloatT fRelTol = 1e-9;
  oops::TestReference::FloatT fAbsTol = 0;

  // Extra parentheses protect EXPECT_* macros arguments which cannot contain commas.
  EXPECT_NO_THROW((oops::TestReference::compare(test1, good_ref1, fRelTol, fAbsTol, iTol)));

  EXPECT_NO_THROW((oops::TestReference::compare(test2, good_ref2, fRelTol, fAbsTol, iTol)));
  EXPECT_NO_THROW((oops::TestReference::compare(test2_1, good_ref2, fRelTol, fAbsTol, iTol)));

  EXPECT_THROWS_AS((oops::TestReference::compare(test1, good_ref2, fRelTol, fAbsTol, iTol)),
                    oops::TestReferenceMissingTestLineError);

  EXPECT_THROWS_AS((oops::TestReference::compare(test2, good_ref1, fRelTol, fAbsTol, iTol)),
                    oops::TestReferenceMissingReferenceLineError);
  EXPECT_THROWS_AS((oops::TestReference::compare(test2_1, good_ref1, fRelTol, fAbsTol, iTol)),
                    oops::TestReferenceMissingReferenceLineError);

  EXPECT_THROWS_AS((oops::TestReference::compare(test1, bad_ref1, fRelTol, fAbsTol, iTol)),
                    oops::TestReferenceIntegerMismatchError);

  EXPECT_THROWS_AS((oops::TestReference::compare(test2, bad_ref2_1, fRelTol, fAbsTol, iTol)),
                    oops::TestReferenceFloatMismatchError);

  EXPECT_THROWS_AS((oops::TestReference::compare(test2_1, bad_ref2_1, fRelTol, fAbsTol, iTol)),
                    oops::TestReferenceFloatMismatchError);

  EXPECT_THROWS_AS((oops::TestReference::compare(test2, bad_ref2_2, fRelTol, fAbsTol, iTol)),
                    oops::TestReferenceTextMismatchError);

  EXPECT_THROWS_AS((oops::TestReference::compare(test3, bad_ref3, fRelTol, fAbsTol, iTol)),
                    oops::TestReferenceTextMismatchError);

  EXPECT_NO_THROW((oops::TestReference::compare(test4, good_ref4, fRelTol, fAbsTol, iTol)));

  EXPECT_THROWS_AS((oops::TestReference::compare(test4, bad_ref4_1, fRelTol, fAbsTol, iTol)),
                    oops::TestReferenceFloatMismatchError);

  EXPECT_THROWS_AS((oops::TestReference::compare(test4, bad_ref4_2, fRelTol, fAbsTol, iTol)),
                    oops::TestReferenceIntegerMismatchError);

  EXPECT_THROWS_AS((oops::TestReference::compare(test4, bad_ref4_3, fRelTol, fAbsTol, iTol)),
                    oops::TestReferenceTextMismatchError);

  // Cannot meet the relative tolerance:
  EXPECT_THROWS_AS((oops::TestReference::compare(test5, good_ref5, fRelTol, fAbsTol, iTol)),
                    oops::TestReferenceFloatMismatchError);

  // But can meet the same absolute tolerance:
  fAbsTol = 1e-9;
  EXPECT_NO_THROW((oops::TestReference::compare(test5, good_ref5, fRelTol, fAbsTol, iTol)));

  EXPECT_THROWS_AS((oops::TestReference::compare(test5, bad_ref5, fRelTol, fAbsTol, iTol)),
                    oops::TestReferenceFloatMismatchError);

  EXPECT_THROWS_AS((oops::TestReference::compare(test5, bad_ref5, fRelTol, fAbsTol, iTol)),
                    oops::TestReferenceFloatMismatchError);
}  // CASE("util/TestReference")

class TestReference : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "test::TestReference";}
  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test

#endif  // TEST_UTIL_TESTREFERENCE_H_
