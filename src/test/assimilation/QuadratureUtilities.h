/*
 * (C) Copyright 2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_ASSIMILATION_QUADUTILS_H_
#define TEST_ASSIMILATION_QUADUTILS_H_

#include <string>
#include <vector>

#include "eckit/testing/Test.h"

#include "oops/../test/TestEnvironment.h"
#include "oops/assimilation/QuadratureRules.h"
#include "oops/assimilation/SLCG.h"
#include "oops/base/DiagonalMatrix.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/FloatCompare.h"

#include "test/assimilation/Vector2D.h"
#include "test/assimilation/Vector3D.h"
#include "test/assimilation/Vector4D.h"

namespace test {

  typedef oops::DiagonalMatrix<Vector2D> Matrix2D;
  typedef oops::DiagonalMatrix<Vector3D> Matrix3D;
  typedef oops::DiagonalMatrix<Vector4D> Matrix4D;

  void test_QuadratureUtilities_GaussLegendre()
  {
    /* Test the Gauss-Legendre quadrature utility by
     * integrating 2*sqrt(1 - x^2) over [-1, 1]. The
     * result should equal pi. */
    
    int quadsize = 20;
    std::vector<double> nodes(quadsize, 0);
    std::vector<double> weights(quadsize, 0);

    oops::gaussLegendre(quadsize, nodes, weights);

    double pi_appx    = 0;
    double min_node   = nodes[0];
    double max_node   = nodes[0];
    double min_weight = weights[0];
    double weight_sum = 0;

    for(int q = 0; q < quadsize; q++) {
      pi_appx    += weights[q]*2*sqrt(1 - pow(nodes[q], 2));
      min_node    = fmin(min_node, nodes[q]);
      max_node    = fmax(max_node, nodes[q]);
      min_weight  = fmin(min_weight, weights[q]);
      weight_sum += weights[q];
    }

    EXPECT(oops::is_close_absolute(pi_appx, 4*atan(1), 1.0e-3));
    EXPECT(oops::is_close_absolute(weight_sum, 2.0, 1.0e-5));
    EXPECT(min_node >= -(1 + 1.0e-5));
    EXPECT(max_node <=  (1 + 1.0e-5));
    EXPECT(min_weight >= 0.);
  }

  void test_QuadratureUtilities_SLCG()
  {
    Vector3D diagA(3, 2, 1);
    Matrix3D A(diagA);
    Vector3D b(1, 1, 1);
    std::vector<double> lambda {0, 1, 2};
    std::vector<Vector3D> X;

    Vector3D X0_expected(1.0/3, 0.5, 1);
    Vector3D X1_expected(.25, 1.0/3, .5);
    Vector3D X2_expected(.2, .25, 1.0/3);
    
    const int maxiter                 = 10;
    const double solution_tolerance   = 1e-9;
    const double reduction_tolerance  = 1e-9;
    const double eigenvalue_tolerance = 1e-3;
    std::vector<double> results       = SLCG(X, A, b, lambda, 3, maxiter, reduction_tolerance);

    EXPECT(oops::is_close_absolute(results[0], 0., reduction_tolerance));
    EXPECT(oops::is_close_absolute(results[1], 1., eigenvalue_tolerance));
    EXPECT(oops::is_close_absolute(results[2], 3., eigenvalue_tolerance));

    EXPECT(oops::is_close_absolute(X[0].x(), X0_expected.x(), solution_tolerance));
    EXPECT(oops::is_close_absolute(X[0].y(), X0_expected.y(), solution_tolerance));
    EXPECT(oops::is_close_absolute(X[0].z(), X0_expected.z(), solution_tolerance));

    EXPECT(oops::is_close_absolute(X[1].x(), X1_expected.x(), solution_tolerance));
    EXPECT(oops::is_close_absolute(X[1].y(), X1_expected.y(), solution_tolerance));
    EXPECT(oops::is_close_absolute(X[1].z(), X1_expected.z(), solution_tolerance));

    EXPECT(oops::is_close_absolute(X[2].x(), X2_expected.x(), solution_tolerance));
    EXPECT(oops::is_close_absolute(X[2].y(), X2_expected.y(), solution_tolerance));
    EXPECT(oops::is_close_absolute(X[2].z(), X2_expected.z(), solution_tolerance));
  }

  void test_QuadratureUtilities_assimilation() {
    /* Test the quadrature rules by using SLCG
     * to update a perturbation in a "toy" data
     * assimilation setup.
     * 
     * B = Diagonal([4, 3, 2, 1])
     * H = [1, 0, 0, 0]
     *     [0, 0, 1, 0]
     * R = [1, 0]
     *     [0, 1] */

    // background covariance
    Vector4D diagB(4, 3, 2, 1);
    Matrix4D B(diagB);

    // background in observation space
    Vector2D diagHBHt(4, 2);
    Matrix2D HBHt(diagHBHt);

    // perturbation adjustment matrix
    Vector4D diagA(sqrt(.2), 1, sqrt(1./3), 1);
    Matrix4D A(diagA);

    // prior and posterior perturbations
    Vector4D b_prior(1, 1, 1, 1);
    Vector4D b_post_target(b_prior);
    A.multiply(b_prior, b_post_target);

    // prior perturbation in observation space
    Vector2D Hb_prior(1, 1);

    // setting up parameters for quadrature and Krylov iteration
    int    quadsize = 10;
    int    maxiter  = 2;
    double slcg_tol = 1e-5;

    std::vector<double> nodes, weights;
    std::vector<Vector2D> X;
    oops::gaussLegendre(quadsize, nodes, weights);
    oops::prepareEAKFQuad(quadsize, nodes, weights);

    // computing the terms in the quadrature
    SLCG(X, HBHt, Hb_prior, nodes, quadsize, maxiter, slcg_tol);
    
    // accumulating the terms
    Vector2D tmp2d(0, 0);
    
    for (int q = 0; q < quadsize; q++) {
      tmp2d.axpy(weights[q], X[q]);
    }

    // tmp4d = Ht*tmp2d
    Vector4D tmp4d(tmp2d.x(), 0, tmp2d.y(), 0);
    Vector4D b_post_appx(0, 0, 0, 0);
    B.multiply(tmp4d, b_post_appx);

    b_post_appx *= -1;
    b_post_appx += b_prior;

    Vector4D error(b_post_target);
    error         -= b_post_appx;
    double err_sq  = dot_product(error, error);

    EXPECT(oops::is_close_absolute(err_sq, 0., 1e-3));
  }

  CASE("assimilation/QuadratureUtilities/gaussLegendre") {
    test_QuadratureUtilities_GaussLegendre();
  }

  CASE("assimilation/QuadratureUtilities/SLCG") {
    test_QuadratureUtilities_SLCG();
  }

  CASE("assimilation/QuadratureUtilities/assimilation") {
    test_QuadratureUtilities_assimilation();
  }

  class QuadratureUtilities : public oops::Test {
   private:
    std::string testid() const override {return "test::QuadratureUtilities";}
    void register_tests() const override {}
    void clear() const override {}
  };

}  // namespace test

#endif  // TEST_ASSIMILATION_QUADUTILS_H_
