/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_TRIDIAGSOLVE_H_
#define OOPS_ASSIMILATION_TRIDIAGSOLVE_H_

#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

#include <vector>

namespace oops {

void TriDiagSolve(const std::vector<double> & diag, const std::vector<double> & sub,
                  const std::vector<double> & rhs, std::vector<double> & sol) {
  const double n = rhs.size();
  ASSERT(sub.size() == n-1);
  ASSERT(diag.size() == n);
  sol.resize(n);
  std::vector<double> c(sub);
  std::vector<double> d(rhs);

  // Forward substitution
  c[0] = c[0] / diag[0];
  d[0] = d[0] / diag[0];
  if (n > 2) {
    for (int ii = 1; ii <= n-2; ++ii) {
      double denom = diag[ii]-sub[ii-1]*c[ii-1];
      c[ii] = c[ii]/denom;
      d[ii] = (d[ii]-sub[ii-1]*d[ii-1])/denom;
    }
  }
  d[n-1] = (d[n-1]-sub[n-2]*d[n-2])/(diag[n-1]-sub[n-2]*c[n-2]);
  // Backward substitution
  sol[n-1] = d[n-1];
  for (int ii = n-2; ii >= 0; --ii) {
    sol[ii] = (d[ii]-(c[ii]*sol[ii+1]));
  }
}

//-------------------------------------------------------------------------------------------------
// Computes T and e1*beta0, gives out ss
void blockTriDiagSolve(const std::vector<Eigen::MatrixXd> & alphas,
                       const std::vector<Eigen::MatrixXd> & betas,
                       const Eigen::MatrixXd & beta0, Eigen::MatrixXd & ss,
                       bool & complexValues,
                       const int members) {
  const int iter = alphas.size();
  complexValues = false;
  Eigen::MatrixXd TT = Eigen::MatrixXd::Zero(iter * members, iter * members);

  for (int ii = 0; ii < iter; ++ii) {
    TT.block(ii*members, ii*members, members, members) = alphas[ii];
    if (ii > 0) TT.block((ii-1)*members, ii*members, members, members) = (betas[ii-1]).transpose();
    if (ii < iter - 1) TT.block((ii+1)*members, ii*members, members, members) = (betas[ii]);
  }

  Eigen::MatrixXd e1b0 = Eigen::MatrixXd::Zero(iter * members, members);
  e1b0.block(0, 0, members, members) = beta0;

  ss = TT.llt().solve(e1b0);

  Eigen::VectorXcd eivals = TT.eigenvalues();
  Eigen::MatrixXd eivalsImag = eivals.imag();
  complexValues = (eivalsImag.array() != 0.0).any();
  if (complexValues) {
    throw eckit::BadValue("T matrix has complex values.");
  }
}


}  // namespace oops

#endif  // OOPS_ASSIMILATION_TRIDIAGSOLVE_H_
