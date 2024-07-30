/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_ASSIMILATION_FULLGMRES_H_
#define TEST_ASSIMILATION_FULLGMRES_H_

#include <string>
#include <vector>

#include "eckit/testing/Test.h"

#include "oops/../test/TestEnvironment.h"
#include "oops/assimilation/FullGMRES.h"
#include "oops/base/DiagonalMatrix.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FloatCompare.h"

#include "test/assimilation/Vector3D.h"

namespace test {

  typedef oops::DiagonalMatrix<Vector3D> Matrix3D;

  void test_FullGMRES_FullGMRES()
  {
    // Solve Ax = b using the FullGMRES algorithm

    // Starting guess
    Vector3D x(1, 2, 3);
    const Vector3D b(30, 15, -20);
    // A is a diagonal matrix
    const Vector3D diagA(10, 5, -5);
    const Matrix3D A(diagA);
    // Deliberately incorrect preconditioning in order to trigger
    // multiple iterations of the minimisation
    const Vector3D diagprecond(1, 1, 1);
    const Matrix3D precond(diagprecond);
    const int maxiter = 10;
    const double tolerance = 1e-9;
    // pq and xy are diagnostics for this minimizer
    std::vector<Vector3D> pq;
    std::vector<Vector3D> xy;

    const double normReduction = FullGMRES(x, b, A, precond, maxiter, tolerance, pq, xy);
    EXPECT(oops::is_close_absolute(x.x(), 3.0, 1.0e-9));
    EXPECT(oops::is_close_absolute(x.y(), 3.0, 1.0e-9));
    EXPECT(oops::is_close_absolute(x.z(), 4.0, 1.0e-9));
    EXPECT(oops::is_close_absolute(normReduction, 0.0, 1.0e-9));
  }

  CASE("assimilation/FullGMRES/FullGMRES") {
    test_FullGMRES_FullGMRES();
  }

  class FullGMRES : public oops::Test {
   public:
    using oops::Test::Test;
   private:
    std::string testid() const override {return "test::FullGMRES";}
    void register_tests() const override {}
    void clear() const override {}
  };

}  // namespace test

#endif  // TEST_ASSIMILATION_FULLGMRES_H_
