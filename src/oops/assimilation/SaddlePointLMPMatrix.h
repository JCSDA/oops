/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_SADDLEPOINTLMPMATRIX_H_
#define OOPS_ASSIMILATION_SADDLEPOINTLMPMATRIX_H_

#include <Eigen/Dense>
#include <iostream>
#include <memory>
#include <vector>

#include <boost/noncopyable.hpp>

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFctWeak.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DualVector.h"
#include "oops/assimilation/SaddlePointPrecondMatrix.h"
#include "oops/assimilation/SaddlePointVector.h"

namespace oops {

  /// The preconditioner for the saddle-point minimizer.
  /*!
   *  The preconditioner is obtained by using low-rank updates.
   *  Let us define the matrices:
   *
   *   \f$ R = [ 0  Rp ],     G = [ Zinv*Rp'         0    ]  \f$
   *   \f$     [ Rq  0 ]          [   0         Zinv'*Rq' ]  \f$
   *
   *   \f$ F = I + G*P0inv*R  \f$
   *
   *  where \f$ P0 \f$ is the initial preconditioner and \f$ I \f$ is the
   *  identity matrix of order 2.
   *
   *  Then the preconditioner is updated from
   *
   *  \f$ Pk = P0inv - P0inv*R*Finv*G*P0inv \f$
   *
   *  The solvers represent matrices as objects that implement a "multiply"
   *  method. This class defines objects that apply the saddle-point matrix.
   */

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
class SaddlePointLMPMatrix : private boost::noncopyable {
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef CostFctWeak<MODEL, OBS>         CostFctWeak_;
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef SaddlePointVector<MODEL, OBS>   SPVector_;
  typedef DualVector<MODEL, OBS>          LagVector_;

 public:
  explicit SaddlePointLMPMatrix(const CostFct_ & j);

  void setup(const std::vector<SPVector_> &, const std::vector<SPVector_> &);

  const int & getk() const {return nvec_;}

  void multiply(const SPVector_ &, SPVector_ &) const;

 private:
  void Pinitmultiply(const SPVector_ &, SPVector_ &) const;

  void Gmultiply(const SPVector_ &, Eigen::VectorXd &) const;

  void Rmultiply(Eigen::VectorXd &, SPVector_ &) const;

  void Fmultiply(Eigen::VectorXd &, Eigen::VectorXd &) const;

  const CostFctWeak_ & j_;
  const bool idmodel_;
  std::vector<SPVector_> xyVEC_;
  std::vector<SPVector_> pqVEC_;
  std::vector<LagVector_> RpVEC_;
  std::vector<CtrlInc_>   RqVEC_;
  std::unique_ptr<SPVector_> spvecinit_;
  int nvec_;
  Eigen::MatrixXd ZMat_;
  Eigen::MatrixXd FMat_;
};

// =============================================================================

template<typename MODEL, typename OBS>
SaddlePointLMPMatrix<MODEL, OBS>::SaddlePointLMPMatrix(const CostFct_ & j)
  : j_(dynamic_cast<const CostFctWeak_ &>(j)),
    idmodel_(false)
{}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
  void SaddlePointLMPMatrix<MODEL, OBS>::setup(const std::vector<SPVector_> & xyVEC,
                                          const std::vector<SPVector_> & pqVEC) {
  xyVEC_.clear();
  pqVEC_.clear();
  RpVEC_.clear();
  RqVEC_.clear();
  for (unsigned ii = 0; ii < xyVEC.size(); ++ii) {
    xyVEC_.push_back(xyVEC[ii]);
    pqVEC_.push_back(pqVEC[ii]);
  }

  nvec_ = xyVEC_.size();

  if (nvec_ != 0) {
    spvecinit_.reset(new SPVector_(xyVEC_[0]));
    (*spvecinit_).zero();

    SPVector_ ww(xyVEC_[0]);
    SPVector_ yy(xyVEC_[0]);

    // RpRqVEC = pqVEC - P0 * xyVEC
    for (unsigned jv = 0; jv < nvec_; ++jv) {
      this->Pinitmultiply(xyVEC_[jv], ww);
      yy = pqVEC_[jv];
      yy -= ww;
      RpVEC_.push_back(yy.lambda());
      RqVEC_.push_back(yy.dx());
    }

    // Z = RpVEC'*xVEC
    ZMat_.resize(nvec_, nvec_);
    double calpha = 10;
    for (unsigned jj = 0; jj < nvec_; ++jj) {
      for (unsigned ii = 0; ii < nvec_; ++ii) {
        ZMat_(ii, jj) = calpha * dot_product(RpVEC_[ii], xyVEC_[jj].lambda());
      }
    }

    // Get F Matrix explicitly
    // Set Identity matrix
    Eigen::MatrixXd Ik;
    Ik = Eigen::MatrixXd::Identity(2*nvec_, 2*nvec_);

    // Apply Fmultiply to identity matrix
    FMat_.resize(2*nvec_, 2*nvec_);
    Eigen::VectorXd output(2*nvec_);
    Eigen::VectorXd input(2*nvec_);
    for (int ii = 0; ii < 2*nvec_; ii++) {
      input = Ik.col(ii);
      Fmultiply(input, output);
      FMat_.col(ii) = output;
    }
  }
}

// =============================================================================
template<typename MODEL, typename OBS>
void SaddlePointLMPMatrix<MODEL, OBS>::Pinitmultiply(const SPVector_ & x,
                                                SPVector_ & z) const {
//  P0 = [D   0   L]
//       [0   R   0]
//       [Lt  0   0]
//
//  where M = I.

  CtrlInc_ ww(x.lambda().dx());

// ADJ block
  j_.runADJ(ww, true);
  z.dx() = x.lambda().dx();
  z.dx() -= ww;

// TLM block
  ww = x.dx();
  j_.runTLM(ww, true);
  z.lambda().zero();
  z.lambda().dx() = x.dx();
  z.lambda().dx() -= ww;

// Diagonal block
  DualVector<MODEL, OBS> diag;
  diag.dx(new CtrlInc_(j_.jb()));
  j_.jb().multiplyB(x.lambda().dx(), diag.dx());
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    diag.append(j_.jterm(jj).multiplyCovar(*x.lambda().getv(jj)));
  }

  z.lambda() += diag;

/*
//  *********************************
//  P0 = [D  0  I]
//       [0  R  0]
//       [I  0  0]
  z.dx() = x.lambda().dx();
  DualVector<MODEL, OBS> diag;
  diag.dx(new CtrlInc_(j_.jb()));
  j_.jb().multiplyB(x.lambda().dx(), diag.dx());
  for (unsigned jj = 0; jj < j_.nterms(); ++jj) {
    diag.append(j_.jterm(jj).multiplyCovar(*x.lambda().getv(jj)));
  }
  z.lambda().zero();
  z.lambda().dx() = x.dx();
  z.lambda() += diag;
//  *********************************
*/
}
// =============================================================================
template<typename MODEL, typename OBS>
void SaddlePointLMPMatrix<MODEL, OBS>::multiply(const SPVector_ & x,
                                           SPVector_ & z) const {
// z = P0inv*x - P0inv*R*Finv*G*P0inv*x

  SaddlePointPrecondMatrix<MODEL, OBS> P0inv(j_);
  SPVector_ svw(x);
  SPVector_ svt(x);

  P0inv.multiply(x, z);
  if (nvec_ != 0) {
    Eigen::VectorXd rhs(2*nvec_);
    Eigen::VectorXd finvx(2*nvec_);
    this->Gmultiply(z, rhs);  // rhs = G*z
    finvx = FMat_.colPivHouseholderQr().solve(rhs);  // finvx = Finv*rhs
    this->Rmultiply(finvx, svw);  // svw = R * finvx
    P0inv.multiply(svw, svt);     // svt = P0inv * svw
    z -= svt;  // z = z - svt
  }
}

// =============================================================================
template<typename MODEL, typename OBS>
void SaddlePointLMPMatrix<MODEL, OBS>::Gmultiply(const SPVector_ & x,
                                            Eigen::VectorXd & z) const {
  Eigen::VectorXd v1(nvec_);
  Eigen::VectorXd v2(nvec_);
  Eigen::VectorXd ax(nvec_);
  Eigen::VectorXd axt(nvec_);
  Eigen::MatrixXd ZMatt(nvec_, nvec_);
  ZMatt = ZMat_.transpose();

  // Rp'*x and Rq'*y
  for (unsigned jj = 0; jj < nvec_; ++jj) {
    v1(jj) = dot_product(RpVEC_[jj], x.lambda());
    v2(jj) = dot_product(RqVEC_[jj], x.dx());
  }

  ax = ZMat_.colPivHouseholderQr().solve(v1);   // Zinv*Rp'*x
  axt = ZMatt.colPivHouseholderQr().solve(v2);  // Ztinv*Rq'*y
  z << ax, axt;  // z = [ax ; axt]
}

// =============================================================================
template<typename MODEL, typename OBS>
void SaddlePointLMPMatrix<MODEL, OBS>::Rmultiply(Eigen::VectorXd & x, SPVector_ & rx) const {
  LagVector_ ww(rx.lambda());
  ww.zero();
  CtrlInc_ zz(rx.dx());
  zz.zero();

  for (unsigned ii = 0; ii < nvec_; ++ii) {
    ww.axpy(x(ii + nvec_), RpVEC_[ii]);
    zz.axpy(x(ii), RqVEC_[ii]);
  }
  rx.lambda() = ww;
  rx.dx() = zz;
}

// =============================================================================
template<typename MODEL, typename OBS>
void SaddlePointLMPMatrix<MODEL, OBS>::Fmultiply(Eigen::VectorXd & x,
                                            Eigen::VectorXd & fx) const {
  SaddlePointPrecondMatrix<MODEL, OBS> P0inv(j_);
  SPVector_ tmp(*spvecinit_);
  SPVector_ tmp2(*spvecinit_);
  Eigen::VectorXd ww(2*nvec_);

  fx = x;
  this->Rmultiply(x, tmp);
  P0inv.multiply(tmp, tmp2);
  this->Gmultiply(tmp2, ww);
  fx += ww;
}

// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_ASSIMILATION_SADDLEPOINTLMPMATRIX_H_
