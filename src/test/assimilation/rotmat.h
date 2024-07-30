/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_ASSIMILATION_ROTMAT_H_
#define TEST_ASSIMILATION_ROTMAT_H_

#include <string>
#include <vector>

#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/assimilation/rotmat.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FloatCompare.h"

namespace test {

  void test_b_zero()
  {
    const double a = 1.0;
    const double b = 0.0;
    double c = -1;
    double s = -1;

    oops::rotmat(a, b, c, s);

    EXPECT_EQUAL(c, 1.0);
    EXPECT_EQUAL(s, 0.0);
  }

  void test_a_zero()
  {
    const double a = 0.0;
    const double b = 1.0;
    double c = -1;
    double s = -1;

    oops::rotmat(a, b, c, s);

    EXPECT_EQUAL(c, 0.0);
    EXPECT_EQUAL(s, 1.0);
  }

  void test_b_gt_a()
  {
    const double a = 3.0;
    const double b = 4.0;
    double c = -1;
    double s = -1;

    oops::rotmat(a, b, c, s);

    EXPECT(oops::is_close_absolute(c, 0.6, 1e-9));
    EXPECT(oops::is_close_absolute(s, 0.8, 1e-9));
  }

  void test_b_leq_a()
  {
    const double a = 4.0;
    const double b = 3.0;
    double c = -1;
    double s = -1;

    oops::rotmat(a, b, c, s);

    EXPECT(oops::is_close_absolute(c, 0.8, 1e-9));
    EXPECT(oops::is_close_absolute(s, 0.6, 1e-9));
  }

  CASE("assimilation/rotmat/b_zero") {
    test_b_zero();
  }

  CASE("assimilation/rotmat/a_zero") {
    test_a_zero();
  }

  CASE("assimilation/rotmat/b_gt_a") {
    test_b_gt_a();
  }

  CASE("assimilation/rotmat/b_leq_a") {
    test_b_leq_a();
  }

class rotmat : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "test::rotmat";}
  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test

#endif  // TEST_ASSIMILATION_ROTMAT_H_
