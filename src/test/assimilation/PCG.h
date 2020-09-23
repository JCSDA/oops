/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_ASSIMILATION_PCG_H_
#define TEST_ASSIMILATION_PCG_H_

#include <string>
#include <vector>

#include "eckit/testing/Test.h"

#include "oops/../test/TestEnvironment.h"
#include "oops/assimilation/IPCG.h"
#include "oops/assimilation/PCG.h"
#include "oops/base/DiagonalMatrix.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FloatCompare.h"

#include "test/assimilation/Vector3D.h"

namespace test {

  void test_Vector3D()
  {
    // Test each function in the Vector3D class

    const Vector3D v1(1.0, 2.0, 3.0);
    const Vector3D v2(4.0, 5.0, 6.0);

    Vector3D vEq = v1;
    EXPECT(vEq.x() == 1.0 && vEq.y() == 2.0 && vEq.z() == 3.0);

    const Vector3D vCopy(v1);
    EXPECT(vCopy.x() == 1.0 && vCopy.y() == 2.0 && vCopy.z() == 3.0);

    Vector3D vAdd = v1;
    vAdd += v2;
    EXPECT(vAdd.x() == 5.0 && vAdd.y() == 7.0 && vAdd.z() == 9.0);

    Vector3D vSub = v2;
    vSub -= v1;
    EXPECT(vSub.x() == 3.0 && vSub.y() == 3.0 && vSub.z() == 3.0);

    Vector3D vMult = v1;
    vMult *= 2.0;
    EXPECT(vMult.x() == 2.0 && vMult.y() == 4.0 && vMult.z() == 6.0);

    Vector3D vMultV = v1;
    vMultV *= v2;
    EXPECT(vMultV.x() == 4.0 && vMultV.y() == 10.0 && vMultV.z() == 18.0);

    Vector3D vDivV = v2;
    vDivV /= v1;
    EXPECT(vDivV.x() == 4.0 && vDivV.y() == 2.5 && vDivV.z() == 2.0);

    Vector3D vAxpy = v1;
    vAxpy.axpy(3, v2);
    EXPECT(vAxpy.x() == 13.0 && vAxpy.y() == 17.0 && vAxpy.z() == 21.0);

    const double dotprod = v1.dot_product_with(v2);
    EXPECT(dotprod == 32.0);
  }

  void test_PCG()
  {
    // Solve Ax = b using the PCG algorithm

    // Starting guess
    Vector3D x(1, 2, 3);
    const Vector3D b(30, 15, -20);
    // A is a diagonal matrix
    const Vector3D diagA(10, 5, -5);
    const oops::DiagonalMatrix<Vector3D> A(diagA);
    // Deliberately incorrect preconditioning in order to trigger
    // multiple iterations of the minimisation
    const Vector3D diagprecond(1, 1, 1);
    const oops::DiagonalMatrix<Vector3D> precond(diagprecond);
    const int maxiter = 10;
    const double tolerance = 1e-9;

    double normReduction = PCG(x, b, A, precond, maxiter, tolerance);
    EXPECT(oops::is_close_absolute(x.x(), 3.0, 1.0e-9));
    EXPECT(oops::is_close_absolute(x.y(), 3.0, 1.0e-9));
    EXPECT(oops::is_close_absolute(x.z(), 4.0, 1.0e-9));
    EXPECT(oops::is_close_absolute(normReduction, 0.0, 1.0e-9));
  }

  void test_IPCG()
  {
    // Solve Ax = b using the IPCG algorithm

    // Starting guess
    Vector3D x(1, 2, 3);
    const Vector3D b(30, 15, -20);
    // A is a diagonal matrix
    const Vector3D diagA(10, 5, -5);
    const oops::DiagonalMatrix<Vector3D> A(diagA);
    // Deliberately incorrect preconditioning in order to trigger
    // multiple iterations of the minimisation
    const Vector3D diagprecond(1, 1, 1);
    const oops::DiagonalMatrix<Vector3D> precond(diagprecond);
    const int maxiter = 10;
    const double tolerance = 1e-9;

    double normReduction = IPCG(x, b, A, precond, maxiter, tolerance);
    EXPECT(oops::is_close_absolute(x.x(), 3.0, 1.0e-9));
    EXPECT(oops::is_close_absolute(x.y(), 3.0, 1.0e-9));
    EXPECT(oops::is_close_absolute(x.z(), 4.0, 1.0e-9));
    EXPECT(oops::is_close_absolute(normReduction, 0.0, 1.0e-9));
  }

  CASE("assimilation/PCG/Vector3D") {
    test_Vector3D();
  }

  CASE("assimilation/PCG/PCG") {
    test_PCG();
  }

  CASE("assimilation/PCG/IPCG") {
    test_IPCG();
  }

  class PCG : public oops::Test {
   private:
    std::string testid() const override {return "test::PCG";}
    void register_tests() const override {}
    void clear() const override {}
  };

}  // namespace test

#endif  // TEST_ASSIMILATION_PCG_H_
