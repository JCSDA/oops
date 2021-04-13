/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_DRPBLOCKLANCZOSMINIMIZER_H_
#define OOPS_ASSIMILATION_DRPBLOCKLANCZOSMINIMIZER_H_

#include <Eigen/Dense>

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

#include "oops/assimilation/BMatrix.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/DRMinimizer.h"
#include "oops/assimilation/HtRinvHMatrix.h"
#include "oops/assimilation/MinimizerUtils.h"
#include "oops/assimilation/TriDiagSolve.h"
#include "oops/mpi/mpi.h"
#include "oops/util/dot_product.h"
#include "oops/util/formats.h"
#include "oops/util/Logger.h"

namespace oops {

template<typename MODEL, typename OBS> class DRPBlockLanczosMinimizer :
         public DRMinimizer<MODEL, OBS> {
  typedef BMatrix<MODEL, OBS>              Bmat_;
  typedef CostFunction<MODEL, OBS>         CostFct_;
  typedef ControlIncrement<MODEL, OBS>     CtrlInc_;
  typedef HtRinvHMatrix<MODEL, OBS>        HtRinvH_;
  typedef Eigen::VectorXd                  eigenvec_;
  typedef Eigen::MatrixXd                  eigenmat_;

 public:
  const std::string classname() const override {return "DRPBlockLanczosMinimizer";}
  DRPBlockLanczosMinimizer(const eckit::Configuration &, const CostFct_ &);
  ~DRPBlockLanczosMinimizer() {}

 private:
  double solve(CtrlInc_ &, CtrlInc_ &, CtrlInc_ &, const Bmat_ &, const HtRinvH_ &, const double,
               const double, const int, const double) override;

  // other functions:
  void get_proj(const CtrlInc_ &, const CtrlInc_ &, eigenmat_ &, int &,
                const eckit::mpi::Comm &, CtrlInc_ &);
  void apply_proj(CtrlInc_ &, const CtrlInc_ &, const eigenmat_ &, int &,
                  const eckit::mpi::Comm &, CtrlInc_ &);
  void mqrgs(CtrlInc_ &, CtrlInc_ &, eigenmat_ &, const CtrlInc_ &, int &,
             const eckit::mpi::Comm &, CtrlInc_ &, CtrlInc_ &);
  void HtRinvH0(const CtrlInc_ &, CtrlInc_ &, const HtRinvH_ &, int &,
                const eckit::mpi::Comm &, CtrlInc_ &);

  // For MPI purposes
  const int members_;
  const int ntasks_;
  const int tasks_per_member_;
  const int global_task_;
  const int mymember_;
  const int local_task_;

  // For diagnostics
  eckit::LocalConfiguration diagConf_;
  int outerLoop_;
};

// ===============================================================================================
// Block-B-Lanczos algorithm
// Primal space
// Lanczos algorithm with complete reorthogonalization.
// MPI version (storage of partial Krylov base on each processor)

template<typename MODEL, typename OBS>
DRPBlockLanczosMinimizer<MODEL, OBS>::DRPBlockLanczosMinimizer(const eckit::Configuration & conf,
                                                      const CostFct_ & J)
  : DRMinimizer<MODEL, OBS>(J), members_(conf.getInt("members")),
    ntasks_(oops::mpi::world().size()),
    tasks_per_member_(ntasks_/members_), global_task_(oops::mpi::world().rank()),
    mymember_(global_task_ / tasks_per_member_), local_task_(global_task_%tasks_per_member_),
    diagConf_(conf), outerLoop_(0) {}

// -----------------------------------------------------------------------------------------------

template<typename MODEL, typename OBS>
double DRPBlockLanczosMinimizer<MODEL, OBS>::solve(CtrlInc_ & xx, CtrlInc_ & xh, CtrlInc_ & rr,
                                            const Bmat_ & B, const HtRinvH_ & HtRinvH,
                                            const double costJ0Jb, const double costJ0JoJc,
                                            const int maxiter, const double tolerance) {
  eigenvec_ zerov = Eigen::VectorXd::Zero(members_);
  eigenmat_ zeromm = Eigen::MatrixXd::Zero(members_, members_);

  CtrlInc_ ww(xh);  // Current w_i
  CtrlInc_ vv(rr);  // Current v_i
  CtrlInc_ zz(xh);  // Current z_i

  CtrlInc_ temp1(xh);  // Current w_i
  CtrlInc_ temp2(xh);  // Current w_i

  int gestag = 0;

  std::vector<std::unique_ptr<CtrlInc_>> Zbase;  // store the zz during iterative process
  std::vector<std::unique_ptr<CtrlInc_>> Vbase;  // store the vv during iterative process

  eigenmat_ projsol(members_, members_);  // projectors (w,z) for residuals calculation.

  eigenmat_ alpha(members_, members_);
  eigenmat_ beta = zeromm;
  eigenmat_ beta0 = zeromm;
  eigenmat_ SSLK;

  std::vector<eigenmat_> ALPHAS;  // store the diagonal blocks of the Arnoldi matrix
  std::vector<eigenmat_> BETAS;  // store the underdiagonal blocks of the Arnoldi matrix

  eigenmat_ ss;  // contains the solution
  eigenvec_ ss_loc = zerov;

  double norm_iiter = 0;
  double norm0 = 0;
  double normReduction = 1;
  double normReductionIter = 1;

  int iterTotal = maxiter;

  bool complexValues = false;

  eigenvec_ norm_red_loc(maxiter);
  eigenmat_ norm_red_all(maxiter, members_);

  eigenvec_ costj_loc(maxiter);
  eigenmat_ costj_all(maxiter, members_);

// -----------------------------------------------------------------------------------------------
// Creating the proper communicator (geographic area) for send and receive
  std::string CommGeoStr = "comm_geo_" + std::to_string(local_task_);
  char const *CommGeoName = CommGeoStr.c_str();
  const eckit::mpi::Comm & CommGeo = oops::mpi::world().split(local_task_, CommGeoName);
  Log::info() << "Geo communicators created, for task 0 of member 0: "
              << CommGeo.name() << " of size " << CommGeo.size() << std::endl;

// -----------------------------------------------------------------------------------------------
  B.multiply(rr, zz);  // z_0 = B * r_0
  norm0 = sqrt(dot_product(zz, rr));

  // QR decomposition
  mqrgs(zz, vv, beta0, rr, gestag, CommGeo, temp1, temp2);  // [z_1, v_1, b0] = qr[r_0, v_0]

  Vbase.emplace_back(std::unique_ptr<CtrlInc_>(new CtrlInc_(vv)));
  Zbase.emplace_back(std::unique_ptr<CtrlInc_>(new CtrlInc_(zz)));

  for (int iiter = 0; iiter < maxiter && normReductionIter > tolerance; ++iiter) {
    Log::info() << "BlockBLanczos starting iteration " << iiter+1 << " for rank: " << mymember_
                << std::endl;

    // Hessian application: w_i = v_i + HtRinvH * B*v_i = v_i + HtRinvH * z_i
    // --> new search directions
    // HtRinvH.multiply(zz, ww);
    HtRinvH0(zz, ww, HtRinvH, gestag, CommGeo, temp1);

    ww += vv;

    // Orthogonalize ww against previous base vectors
    for (int jiter = 0; jiter < iiter + 1; ++jiter) {
      get_proj(ww, *Zbase[jiter], alpha, gestag, CommGeo, temp1);
      apply_proj(ww, *Vbase[jiter], alpha, gestag, CommGeo, temp1);
    }

    B.multiply(ww, zz);
    get_proj(ww, zz, projsol, gestag, CommGeo, temp1);  // projectors for residuals calculation

    vv = ww;
    mqrgs(zz, vv, beta, ww, gestag, CommGeo, temp1, temp2);

    Zbase.emplace_back(std::unique_ptr<CtrlInc_>(new CtrlInc_(zz)));
    Vbase.emplace_back(std::unique_ptr<CtrlInc_>(new CtrlInc_(vv)));
    ALPHAS.push_back(alpha);
    BETAS.push_back(beta);

    // solve T ss = beta0 * e1
    blockTriDiagSolve(ALPHAS, BETAS, beta0, ss, complexValues, members_);

    eigenvec_ ss_loc = (ss.block(iiter*members_, 0, members_, members_)).col(mymember_);

    eigenvec_ temp = projsol * ss_loc;
    norm_iiter = sqrt(temp.dot(ss_loc));
    normReduction = norm_iiter / norm0;

    const double costJ0 = costJ0Jb + costJ0JoJc;
    double costJ = costJ0;
    double costJb = 0;
    double costJoJc = 0;

    for (int ll = 0; ll < iiter+1; ++ll) {
      SSLK = - (ss.block(ll*members_, 0, members_, members_));
      // Compute the quadratic cost function
      // J[du_{i}] = J[0] - 0.5 s_{i}^T Z_{i}^T r_{0}
      // Jb[du_{i}] = 0.5 s_{i}^T V_{i}^T Z_{i} s_{i}
      temp2.zero();
      apply_proj(temp2, *Zbase[ll], SSLK, gestag, CommGeo, temp1);
      costJ -= 0.5 * dot_product(temp2, rr);
    }

    Log::info() << "BlockBLanczos end of iteration " << iiter+1 << std::endl;
    printNormReduction(iiter+1, norm_iiter, normReduction);
    printQuadraticCostFunction(iiter+1, costJ, costJb, costJoJc);

    norm_red_loc[iiter] = normReduction;
    costj_loc[iiter] = costJ;

    oops::mpi::world().allReduce(normReduction, normReductionIter, eckit::mpi::max());
    if (normReductionIter < tolerance) {
      Log::info() << "DRPBlockLanczos: Achieved required reduction in residual norm." << std::endl;
      iterTotal = iiter+1;
    }
  }  // main loop iiter

  oops::mpi::allGather(CommGeo, norm_red_loc, norm_red_all);
  oops::mpi::allGather(CommGeo, costj_loc, costj_all);

  xh.zero();
  xx.zero();
  eigenmat_ tmp_norm;
  eigenmat_ tmp_costj;
  Eigen::IOFormat HeavyFmt(-2, 0, ", ", ";\n");
  Eigen::IOFormat TestFmt(-1, 0, ", ", ";\n");


  for (int ll = 0; ll < iterTotal; ++ll) {
    SSLK = - (ss.block(ll*members_, 0, members_, members_));
    apply_proj(xh, *Vbase[ll], SSLK, gestag, CommGeo, temp1);
    apply_proj(xx, *Zbase[ll], SSLK, gestag, CommGeo, temp1);
    if (outerLoop_ == 0) writeKrylovBasis(diagConf_, *Zbase[ll], ll);
    tmp_norm = norm_red_all.block(ll, 0, 1, members_);
    tmp_costj = costj_all.block(ll, 0, 1, members_);

    Log::info() << "   Norm reduction all members (" << std::setw(2) << ll+1 << ") = "
                << tmp_norm.format(HeavyFmt) << std::endl;
    Log::info() << "   Quadratic cost function all members: J (" << std::setw(2) << ll+1 << ") = "
                << tmp_costj.format(HeavyFmt) << std::endl;
    Log::test() << "   Norm reduction all members (" << std::setw(2) << ll+1 << ") = "
                << tmp_norm.format(TestFmt) << std::endl;
    Log::test() << "   Quadratic cost function all members: J (" << std::setw(2) << ll+1 << ") = "
                << tmp_costj.format(TestFmt) << std::endl;
  }

  eckit::mpi::deleteComm(CommGeoName);
  ++outerLoop_;
  return normReduction;
}

// -----------------------------------------------------------------------------------------------
//                                       OTHER FUNCTIONS
// -----------------------------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void DRPBlockLanczosMinimizer<MODEL, OBS>::get_proj(const CtrlInc_ & www,
                                             const CtrlInc_ & incr_tosend,
                                             eigenmat_ & alpha_mat, int & tag,
                                             const eckit::mpi::Comm & comm,
                                             CtrlInc_ & incr_rcv) {
  // Computes mat_proj = Zt.W
  // with incr_tosend(member_i) the columns of matrix Z
  // with www(member_i) the columns of matrix W
  eigenvec_ alpha_loc = Eigen::VectorXd::Zero(members_);
  int tag_rcv;
  int tag_send;

  for (int p = 0; p < mymember_; p++) {
    tag_rcv = p * members_ + mymember_ + tag * members_ * members_;
    tag_send = mymember_ * members_ + p + tag * members_ * members_;
    oops::mpi::send(comm, incr_tosend, p, tag_send);
    oops::mpi::receive(comm, incr_rcv, p, tag_rcv);
    alpha_loc(p) = dot_product(incr_rcv, www);
  }

  alpha_loc(mymember_) = dot_product(incr_tosend, www);

  for (int p = mymember_ + 1; p < members_; p++) {
    tag_rcv = p * members_ + mymember_ + tag * members_ * members_;
    tag_send = mymember_ * members_ + p + tag * members_ * members_;
    oops::mpi::receive(comm, incr_rcv, p, tag_rcv);
    oops::mpi::send(comm, incr_tosend, p, tag_send);
    alpha_loc(p) = dot_product(incr_rcv, www);
  }
  oops::mpi::allGather(comm, alpha_loc, alpha_mat);

  tag += 2;
}

// -----------------------------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void DRPBlockLanczosMinimizer<MODEL, OBS>::apply_proj(CtrlInc_ & incr_tochange,
                                               const CtrlInc_ & incr_tosend,
                                               const eigenmat_ & alpha_mat, int & tag,
                                               const eckit::mpi::Comm & comm,
                                               CtrlInc_ & incr_rcv) {
  // Computes W = W - V.alpha
  // with incr_tochange(member_i) the columns of matrix W
  // with incr_tosend(member_i) the columns of matrix V
  eigenvec_ alpha_loc = alpha_mat.col(mymember_);
  int tag_rcv;
  int tag_send;

  for (int p = 0; p < mymember_; p++) {
    tag_rcv = p * members_ + mymember_ + tag * members_ * members_;
    tag_send = mymember_ * members_ + p + tag * members_ * members_;
    oops::mpi::send(comm, incr_tosend, p, tag_send);
    oops::mpi::receive(comm, incr_rcv, p, tag_rcv);
    incr_tochange.axpy(-alpha_loc(p), incr_rcv);
  }

  incr_tochange.axpy(-alpha_loc(mymember_), incr_tosend);

  for (int p = mymember_ + 1; p < members_; p++) {
    tag_rcv = p * members_ + mymember_ + tag * members_ * members_;
    tag_send = mymember_ * members_ + p + tag * members_ * members_;
    oops::mpi::receive(comm, incr_rcv, p, tag_rcv);
    oops::mpi::send(comm, incr_tosend, p, tag_send);
    incr_tochange.axpy(-alpha_loc(p), incr_rcv);
  }
  tag += 2;
}

// -----------------------------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void DRPBlockLanczosMinimizer<MODEL, OBS>::mqrgs(CtrlInc_ & zzz, CtrlInc_ & vvv,
                                         eigenmat_ & beta_mat, const CtrlInc_ & www,
                                         int & tag,
                                         const eckit::mpi::Comm & comm,
                                         CtrlInc_ & v_other, CtrlInc_ & z_other) {
  // QR decomposition: [zzz, vvv, beta_mat] = qr[www, vvv] using Gram-Schmidt
  eigenvec_ beta_loc = Eigen::VectorXd::Zero(members_);
  int tag_rcv_v;
  int tag_rcv_z;
  int tag_send_v;
  int tag_send_z;

  // Orthogonalization
  for (int p = 0; p < mymember_; p++) {
    tag_rcv_v = p * members_ + mymember_ + tag * members_ * members_;
    tag_rcv_z = p * members_ + mymember_ + (tag + 1) * members_ * members_;
    oops::mpi::receive(comm, v_other, p, tag_rcv_v);
    oops::mpi::receive(comm, z_other, p, tag_rcv_z);
    beta_loc(p) = dot_product(www, z_other) / dot_product(v_other, z_other);
    vvv.axpy(-beta_loc(p), v_other);
    zzz.axpy(-beta_loc(p), z_other);
  }

  for (int p = mymember_ + 1; p < members_; p++) {
    tag_send_v = mymember_ * members_ + p + tag * members_ * members_;
    tag_send_z = mymember_ * members_ + p + (tag + 1) * members_ * members_;
    oops::mpi::send(comm, vvv, p, tag_send_v);
    oops::mpi::send(comm, zzz, p, tag_send_z);
  }

  // Normalization (TODO: needs better writing)
  beta_loc(mymember_) = sqrt(std::max(dot_product(vvv, zzz), 1e-15));

  vvv *= (1 / beta_loc(mymember_));
  zzz *= (1 / beta_loc(mymember_));

  oops::mpi::allGather(comm, beta_loc, beta_mat);

  for (int ii = 0; ii < members_; ++ii) {
    for (int jj = ii + 1; jj < members_; ++jj) {
      beta_mat(ii, jj) =  beta_mat(ii, jj) * beta_mat(ii, ii);
    }
  }

  tag += 2;
}

// -----------------------------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void DRPBlockLanczosMinimizer<MODEL, OBS>::HtRinvH0(const CtrlInc_ & z_loc, CtrlInc_ & w_out,
                                             const HtRinvH_ & HtRinvH, int & tag,
                                             const eckit::mpi::Comm & comm, CtrlInc_ & z_other) {
// send z_loc to task 0, process HtRinvH_0 * z_loc and send the result back to original task
  int tag_send_w;
  int tag_rcv_w;
  int tag_send_z;
  int tag_rcv_z;
  int dest = 0;

  if (mymember_ == 0) {
    for (int ii = 1; ii < members_; ++ii) {
      tag_rcv_z = members_ * ii + tag * members_ * members_;
      oops::mpi::receive(comm, z_other, ii, tag_rcv_z);
      HtRinvH.multiply(z_other, w_out);
      tag_send_w = members_ * ii + (tag + 1) * members_ * members_;
      oops::mpi::send(comm, w_out, ii, tag_send_w);  // send to task ii
    }
    HtRinvH.multiply(z_loc, w_out);
  } else {
      tag_send_z = mymember_ * members_ + tag * members_ * members_;
      oops::mpi::send(comm, z_loc, dest, tag_send_z);  // send to task 0
      tag_rcv_w = mymember_ * members_ + (tag + 1) * members_ * members_;
      oops::mpi::receive(comm, w_out, dest, tag_rcv_w);
  }
  tag += 2;
}

// -----------------------------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_DRPBLOCKLANCZOSMINIMIZER_H_
