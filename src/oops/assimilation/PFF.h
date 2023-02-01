/*
 * (C) Copyright 2022-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_PFF_H_
#define OOPS_ASSIMILATION_PFF_H_

#include <iomanip>
#include <memory>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/mpi/Comm.h"

#include "oops/assimilation/BMatrix.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostFunction.h"
#include "oops/assimilation/Minimizer.h"
#include "oops/base/State.h"

#include "oops/util/Logger.h"

namespace oops {

/// PFF
/*!
 * Implements the Particle Flow Filter using MPI.
 * A particle flow filter for high-dimensional system applications,
 * Chih-Chi Hu, Peter Jan Van Leeuwen, QJRMS 2021
 *
 * standard deviation is set in the yaml file
 */

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS> class PFF : public Minimizer<MODEL, OBS> {
  typedef BMatrix<MODEL, OBS>             Bmat_;
  typedef ControlIncrement<MODEL, OBS>    CtrlInc_;
  typedef ControlVariable<MODEL, OBS>     CtrlVar_;
  typedef State<MODEL>                    State_;
  typedef CostFunction<MODEL, OBS>        CostFct_;
  typedef Minimizer<MODEL, OBS>           Minimizer_;

 public:
  const std::string classname() const override {return "PFF";}
  PFF(const eckit::Configuration & conf, const CostFct_ & J)
    : Minimizer_(J), J_(J), gradJb_(), minimConf_(conf) {}
  ~PFF() {}

 private:
  CtrlInc_ * doMinimize(const eckit::Configuration &) override;
  CtrlInc_  kernel(const CtrlVar_ &, const CtrlVar_ &);
  CtrlInc_  doExponential(const CtrlInc_ &);
  void computeMean(const CtrlVar_ &, CtrlVar_ &, CtrlVar_ &, const eckit::mpi::Comm &);
  void computeNorm(CtrlInc_ &, CtrlInc_ &, const eckit::mpi::Comm &);
  void diagonalsOfB(const Bmat_ &, double);
  void dxi(CtrlInc_ &, const CtrlVar_ &, CtrlVar_ &, CtrlVar_ &, const Bmat_ &,
          CtrlInc_ &, CtrlInc_ &, const eckit::mpi::Comm &);
  double norm(CtrlInc_);
  int CtrlIncToVecIndex(int);
  int generateTag(int, int);

  const CostFct_ & J_;
  std::unique_ptr<CtrlInc_> gradJb_;
  eckit::LocalConfiguration minimConf_;
  double diagB_;
  double invB_;
  std::vector<double> norms = {1, 1, 1};

  int iter_;
  int ct_ = 0;
  int ctCheck_;
  int myrank_;
  int nmembers_;
  int resolution_;
  double norm_;
  double weight_;
  double eps_;
  double minLearningRate_;
  double inflationFactor_;
};

// =============================================================================

template<typename MODEL, typename OBS>
ControlIncrement<MODEL, OBS> *
PFF<MODEL, OBS>::doMinimize(const eckit::Configuration & config) {
  iter_ = config.getInt("iteration");
  resolution_ = config.getInt("geometry.resol");
  eps_ = minimConf_.getDouble("eps");
  minLearningRate_ = minimConf_.getDouble("minimum learning rate");
  ctCheck_ = minimConf_.getDouble("ct check");
  inflationFactor_ = minimConf_.getDouble("inflation factor");
  int ninner = config.getInt("ninner");
  double gnreduc = config.getDouble("gradient norm reduction");
  double standardDev = minimConf_.getDouble("standard_deviation");

  if (gradJb_) {
    gradJb_.reset(new CtrlInc_(J_.jb().resolution(), *gradJb_));
  } else {
    gradJb_.reset(new CtrlInc_(J_.jb()));
  }

  Log::info() << "PFF: max iter = " << ninner
              << ", requested norm reduction = " << gnreduc << std::endl;


// Define the matrices:
  const Bmat_ B(J_);

// Define minimisation starting point
  // dx
  CtrlInc_ * dx = new CtrlInc_(J_.jb());
  CtrlInc_ * rr = new CtrlInc_(J_.jb().getFirstGuess());
  CtrlInc_ rr__(J_.jb().getFirstGuess());  // remove later
  CtrlInc_ rrOther(rr__);
  CtrlInc_ zz(J_.jb());

  J_.computeGradientFG(*rr);
  J_.jb().addGradientFG(*rr, *gradJb_);

// Communicator and MPI parameters:
  const eckit::mpi::Comm & comm = oops::mpi::world();
  myrank_ = comm.rank();
  nmembers_ = comm.size();
  weight_ = 1/static_cast<double>(nmembers_);
  Log::info() << "Running the PFF algorithm with " << nmembers_ << " particles" << std::endl;

// Other useful terms
  CtrlVar_ xx(J_.jb().getBackground());  // get initial state
  CtrlVar_ xxOther(xx);  // to use when receiving
  CtrlVar_ xxMean(xx);

  computeMean(xx, xxOther, xxMean, comm);
  diagonalsOfB(B, standardDev);
  dxi(*dx, xx, xxOther, xxMean, B, *rr, rrOther, comm);

  Log::info() << classname() << " output increment" << *dx << std::endl;

  return dx;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void PFF<MODEL, OBS>::diagonalsOfB(const Bmat_ & B, double standardDev) {
    double variance = standardDev * standardDev;
    diagB_ = variance;
    invB_ = 1.0/variance;
}
// -----------------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void PFF<MODEL, OBS>::computeMean(const CtrlVar_ & xxlocal, CtrlVar_ & xxother,
                                  CtrlVar_ & xxmean, const eckit::mpi::Comm & comm) {
  // Compute the mean of all states on task 0:
  int tag_send;
  int tag_rcv;
  int dest = 0;

  if (myrank_== 0) {
    State_ xxmeanstate(xxlocal.state());
    xxmeanstate.zero();
    xxmeanstate.accumul(weight_, xxlocal.state());

    for (int ii = 1; ii < nmembers_; ++ii) {
      tag_rcv = 123;
      tag_send = 234;
      oops::mpi::receive(comm, xxother, ii, tag_rcv);
      xxmeanstate.accumul(weight_, xxother.state());
    }

    xxmean.state() = xxmeanstate;

    for (int ii = 1; ii < nmembers_; ++ii) {
      oops::mpi::send(comm, xxmean, ii, tag_send);
    }
  } else {
    tag_rcv = 234;
    tag_send = 123;
    oops::mpi::send(comm, xxlocal, dest, tag_send);
    oops::mpi::receive(comm, xxmean, dest, tag_rcv);
  }
}

// -----------------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void PFF<MODEL, OBS>::dxi(CtrlInc_ & dx, const CtrlVar_ & xxlocal, CtrlVar_ & xxother,
                            CtrlVar_ & xxmean, const Bmat_ & B, CtrlInc_ & rr, CtrlInc_ & rrother,
                            const eckit::mpi::Comm & comm) {
  double alpha = 1/static_cast<double>(nmembers_);

  int tagSend1;
  int tagRcv1;
  int tagSend2;
  int tagRcv2;

  CtrlInc_ rrlocal(rr);
  CtrlInc_ diffState(rrlocal);  // to compute the state difference
  CtrlInc_ kernelNN(rrlocal);
  CtrlInc_ dxij(dx);
  CtrlInc_ gradLogPostK(rr);
  CtrlInc_ xxMeanDiff(rr);

  dxij.zero();

  int i = myrank_;
  rrlocal *= -1.0;
  for (int j = 0; j < nmembers_; j++) {
    gradLogPostK.zero();
    if (myrank_ != j) {
      tagSend1 = generateTag(i, j);
      tagRcv1 = generateTag(j, i);
      oops::mpi::send(comm, xxlocal, j, tagSend1);
      oops::mpi::receive(comm, xxother, j, tagRcv1);

      tagSend2 = generateTag(i, j+1);
      tagRcv2 = generateTag(j, i+1);
      oops::mpi::send(comm, rrlocal, j, tagSend2);
      oops::mpi::receive(comm, rrother, j, tagRcv2);

      diffState.diff(xxlocal, xxother);  // x_i - x_j type CtrlInc_
      kernelNN = kernel(xxlocal, xxother);
      // K(x_i, x_j)o[BH^{T}R^{-1}o[ (y_j - H(x_j)) + (x_j - \mean(x)) ]
      // - B o diag(B^{-1})o(x_i-x_j)oK(x_i, x_j)

      B.multiply(rrother, gradLogPostK);  // B * rr = BH^T R^{-1} (y_j - H(x_j))
      gradLogPostK *= inflationFactor_;
      xxMeanDiff.diff(xxother, xxmean);

      gradLogPostK -= xxMeanDiff;
      gradLogPostK.schur_product_with(kernelNN);

      diffState.schur_product_with(kernelNN);
      diffState *= 1/alpha;
      diffState *= inflationFactor_;
      gradLogPostK -= diffState;

    } else {
      kernelNN = kernel(xxlocal, xxlocal);  // kernel
      // K(x_i, x_j) o B[ H^TR^{-1}(y_j - H(x_j) - (x_j - x_mean)]
      // - B o diag(B)^{-1}o(x_i-x_j)oK(x_i, x_j)
      B.multiply(rrlocal, gradLogPostK);  // BH^TR^{-1}(y_j - H(x_j))
      gradLogPostK *= inflationFactor_;

      xxMeanDiff.diff(xxlocal, xxmean);  // x_j - x_mean

      gradLogPostK -= xxMeanDiff;
      // K(x_i, x_j) o B[ <H^TR^{-1}(y_j-H(x_j))> - <x_j-x_mean> ]
      gradLogPostK.schur_product_with(kernelNN);
    }
    dxij += gradLogPostK;
  }

  dxij *= 1/static_cast<double>(nmembers_);

  CtrlInc_ dxijOther(dxij);

  computeNorm(dxij, dxijOther, comm);
  norms[2] = norm_;

  if (iter_ == 0) {
    norms[0] = norms[2];
    norms[1] = norms[2];
    dxij *= eps_;
    dx = dxij;
    ct_++;
  } else if (eps_ < minLearningRate_) {
    dx.zero();
    Log::info() << "learning rate is too small" << std::endl;
  } else if (iter_ >= 1 && norms[2] > 1.02*norms[1]) {
    eps_ = eps_/1.5;
    iter_--;
    ct_ = 0;
    dx.zero();
    Log::info() << "eps has changed to: " << eps_ << " - redo pff" << std::endl;
  } else if (ct_ >= ctCheck_ && iter_ >= 1 && norms[2] <= 1.02*norms[1]) {
    ct_ = 0;
    dxij *= (eps_ * 1.5);
  } else {
    dxij *= eps_;
    dx = dxij;
    ct_++;
  }
  Log::info() << "iteration: " << iter_ << ", eps: " << eps_ << " norm: "
              << (norms[2]/norms[0]) * 100;
  Log::info() << " norms: [" << norms[0] << ", " << norms[1] << ", "
              << norms[2] << "]" << std::endl;
  Log::test() << "iteration: " << iter_ << ", eps: " << eps_ << " norm: "
              << (norms[2]/norms[0]) * 100;
  Log::test() << " norms: [" << norms[0] << ", " << norms[1] << ", "
              << norms[2] << "]" << std::endl;

  norms[1] = norms[2];

  dx = dxij;
}

// -----------------------------------------------------------------------------------

template<typename MODEL, typename OBS>
double PFF<MODEL, OBS>::norm(CtrlInc_ val) {
  return sqrt(dot_product(val, val));
}

// -----------------------------------------------------------------------------------

template<typename MODEL, typename OBS>
ControlIncrement<MODEL, OBS> PFF<MODEL, OBS>::kernel(const CtrlVar_ & xxLocal,
                                                const CtrlVar_ & xxOther) {
  CtrlInc_ stateDiff(J_.jb());
  CtrlInc_ Knn(stateDiff);
  double alpha = 1/static_cast<double>(nmembers_);

  stateDiff.zero();
  stateDiff.diff(xxLocal, xxOther);  // x_i - x_j
  stateDiff.schur_product_with(stateDiff);  // (x_i - x_j)^2
  stateDiff *= -1/(2.0*alpha*diagB_);  // -1/(2.0*alpha*diagsOfBVec_[CtrlIncToVecIndex(0)]);
  Knn = doExponential(stateDiff);

  return Knn;
}

// -----------------------------------------------------------------------------------

template<typename MODEL, typename OBS>
int PFF<MODEL, OBS>::CtrlIncToVecIndex(int index) {
  int start = 0;
  int vecIndex = index + start;

  return vecIndex;
}

// -----------------------------------------------------------------------------------

template<typename MODEL, typename OBS>
int PFF<MODEL, OBS>::generateTag(int i, int j) {
  int num;
  num = 0.5 * (i + j)*(i + j + 1) + j;

  return num;
}

// -----------------------------------------------------------------------------------

template<typename MODEL, typename OBS>
ControlIncrement<MODEL, OBS>  PFF<MODEL, OBS>::doExponential(const CtrlInc_ & item) {
  CtrlInc_ item_(item);
  std::vector<double> exponentialVec;
  std::vector<double> itemVec;
  item.serialize(itemVec);

  for (int i = 0; i < resolution_; i++) {
    itemVec[i] = exp(itemVec[i]);
  }

  size_t ii = 0;
  item_.deserialize(itemVec, ii);

  return item_;
}

// -----------------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void PFF<MODEL, OBS>::computeNorm(CtrlInc_ & dxij, CtrlInc_ & dxijOther,
                                const eckit::mpi::Comm & comm) {
  std::vector<double> vec;
  int tagSend;
  int tagRcv;
  int i = 0;

  dxijOther.zero();
  if (myrank_ == 0) {
    for (int j = 1; j < nmembers_; j++) {
      // receive
      tagRcv = generateTag(j, i);
      oops::mpi::receive(comm, dxijOther, j, tagRcv);
      dxij += dxijOther;
    }
  } else {
      // send
      tagSend = generateTag(myrank_, i);
      oops::mpi::send(comm, dxij, i, tagSend);
  }

  norm_ = sqrt(dot_product(dxij, dxij) / (static_cast<double>(nmembers_  * resolution_)));
}
// -----------------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_ASSIMILATION_PFF_H_
