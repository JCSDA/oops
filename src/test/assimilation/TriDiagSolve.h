/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_ASSIMILATION_TRIDIAGSOLVE_H_
#define TEST_ASSIMILATION_TRIDIAGSOLVE_H_

#include <Eigen/Dense>

#include <string>
#include <vector>

#include "eckit/testing/Test.h"

#include "oops/../test/TestEnvironment.h"
#include "oops/assimilation/TriDiagSolve.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FloatCompare.h"

namespace test {

  void test_TriDiagSolve()
  {
    const std::vector<double> diag{1.0, 1.0, 1.0};
    const std::vector<double> sub{0.5, 0.5};
    const std::vector<double> rhs{1.0, 0.5, 0.0};
    std::vector<double> sol;

    oops::TriDiagSolve(diag, sub, rhs, sol);

    EXPECT_EQUAL(sol[0], 1.0);
    EXPECT_EQUAL(sol[1], 0.0);
    EXPECT_EQUAL(sol[2], 0.0);
  }

  void test_blockTriDiagSolve()
  {
    const int members = 3;
    Eigen::MatrixXd A0 = Eigen::MatrixXd::Ones(members, members);
    Eigen::MatrixXd A1 = Eigen::MatrixXd::Ones(members, members);
    Eigen::MatrixXd B0 = Eigen::MatrixXd::Ones(members, members);
    Eigen::MatrixXd B1 = Eigen::MatrixXd::Ones(members, members);
    B0(1, 0) = 0;
    B0(2, 0) = 0;
    B0(2, 1) = 0;
    B1(1, 0) = 0;
    B1(2, 0) = 0;
    B1(2, 1) = 0;

    Eigen::MatrixXd beta0 = Eigen::MatrixXd::Ones(members, members);
    beta0(1, 0) = 0;
    beta0(2, 0) = 0;
    beta0(2, 1) = 0;
    Eigen::MatrixXd ss;
    bool complexValues = false;

    oops::blockTriDiagSolve({A0, A1}, {B0, B1}, beta0, ss, complexValues, members);

    EXPECT(oops::is_close_absolute(ss(0, 0), 4.0, 1e-9));
    EXPECT(oops::is_close_absolute(ss(0, 1), 2.0, 1e-9));
    EXPECT(oops::is_close_absolute(ss(0, 2), 2.0, 1e-9));
    EXPECT(oops::is_close_absolute(ss(1, 0), -2.0, 1e-9));
    EXPECT(oops::is_close_absolute(ss(1, 1), 1.0, 1e-9));
    EXPECT(oops::is_close_absolute(ss(1, 2), 0.0, 1e-9));
    EXPECT(oops::is_close_absolute(ss(2, 0), 0.0, 1e-9));
    EXPECT(oops::is_close_absolute(ss(2, 1), -1.0, 1e-9));
    EXPECT(oops::is_close_absolute(ss(2, 2), 1.0, 1e-9));
    EXPECT(oops::is_close_absolute(ss(3, 0), -1.0, 1e-9));
    EXPECT(oops::is_close_absolute(ss(3, 1), -1.0, 1e-9));
    EXPECT(oops::is_close_absolute(ss(3, 2), -2.0, 1e-9));
    EXPECT(oops::is_close_absolute(ss(4, 0), 2.0, 1e-9));
    EXPECT(oops::is_close_absolute(ss(4, 1), 1.0, 1e-9));
    EXPECT(oops::is_close_absolute(ss(4, 2), 1.0, 1e-9));
    EXPECT(oops::is_close_absolute(ss(5, 0), -1.0, 1e-9));
    EXPECT(oops::is_close_absolute(ss(5, 1), 0.0, 1e-9));
    EXPECT(oops::is_close_absolute(ss(5, 2), 0.0, 1e-9));
    EXPECT_EQUAL(complexValues, false);
  }

  CASE("assimilation/TriDiagSolve/TriDiagSolve") {
    test_TriDiagSolve();
  }

  CASE("assimilation/TriDiagSolve/blockTriDiagSolve") {
    test_blockTriDiagSolve();
  }

  class TriDiagSolve : public oops::Test {
   public:
    using oops::Test::Test;
   private:
    std::string testid() const override {return "test::TriDiagSolve";}
    void register_tests() const override {}
    void clear() const override {}
  };

}  // namespace test

#endif  // TEST_ASSIMILATION_TRIDIAGSOLVE_H_
