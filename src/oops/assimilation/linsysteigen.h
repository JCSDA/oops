/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_LINSYSTEIGEN_H_
#define OOPS_ASSIMILATION_LINSYSTEIGEN_H_

#include <Eigen/Dense>

namespace oops {

typedef Eigen::MatrixXd  eigenmat_;
typedef Eigen::VectorXd  eigenvec_;

eigenmat_ linsysteigen(const eigenmat_ & AA, const eigenmat_ & YY) {
  // Solve AX=Y (unknown X), A full square matrix, Y matrix
  return AA.fullPivLu().solve(YY);
}

eigenvec_ linsysteigen(const eigenmat_ & AA, const eigenvec_ & yy) {
  // Solve Ax=y (unknown x), A full square matrix, y vector
  return AA.fullPivLu().solve(yy);
}

eigenmat_ linsysteigentrans(const eigenmat_ & UU, const eigenmat_ & II) {
  // Solve LU=I (unknown L), U & I full square matrixes
  eigenmat_ AA = UU.transpose();
  eigenmat_ BB = II.transpose();
  return linsysteigen(AA, BB).transpose();
}

}       // namespace oops
#endif  // OOPS_ASSIMILATION_LINSYSTEIGEN_H_
