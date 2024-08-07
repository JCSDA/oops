/*
 * (C) Copyright 2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_QUADRATURERULES_H_
#define OOPS_ASSIMILATION_QUADRATURERULES_H_

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>
#include <vector>

#include "oops/util/Logger.h"

namespace oops {

//-------------------------------------------------------------------------------------------------
// Computes nodes and weights for Gauss-Legendre quadrature using the Golub-Welsch algorithm.
void gaussLegendre(const int quadsize, std::vector<double>& nodes, std::vector<double>& weights) {
  Log::info() << "GaussLegendre: Starting and forming Golub-Welsch matrix." << std::endl;

  Eigen::MatrixXf GW(quadsize, quadsize);
  GW.setZero();

  for(int q = 0; q + 1 < quadsize; q++) {
    GW(q+1, q) = .5/sqrt(1 - pow(2*(q + 1), -2));
    GW(q, q+1) = GW(q+1, q);
  }

  Log::info() << "GaussLegendre: Eigenvalue decomposition of Golub-Welsch matrix." << std::endl;

  Eigen::SelfAdjointEigenSolver<Eigen::MatrixXf> solver(quadsize);
  solver.compute(GW);
  Eigen::VectorXf evals = solver.eigenvalues();
  Eigen::MatrixXf evecs = solver.eigenvectors();

  Log::info() << "GaussLegendre: Calculating quadrature weights and nodes." << std::endl;

  nodes.clear();
  weights.clear();

  for(int q = 0; q < quadsize; q++) {
    nodes.push_back(evals(q));
    weights.push_back(2*pow(evecs(0, q), 2));
  }
}

void prepareEAKFQuad(const int quadsize, std::vector<double>& nodes, std::vector<double>& weights) {
  const double PI = 4*atan(1);
  double s, w;

  for (int q = 0; q < quadsize; q++) {
    s = pow(tan(.25*PI*(nodes[q] + 1)), 2);
    w = .5*weights[q]/pow(cos(.25*PI*(nodes[q] + 1)), 2);

    nodes[q]   = s + 1;
    weights[q] = w/(s + 1);
  }
}

}  // namespace oops

#endif  // OOPS_ASSIMILATION_QUADRATURERULES_H_
