/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_QNEWTONLMP_H_
#define OOPS_ASSIMILATION_QNEWTONLMP_H_

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"

namespace oops {

  /*  The QUASI-NEWTON  Limited Memory preconditioner matrix \f$ C \f$
   *  for the system matrix \f$ AB = (I + Ht Rinv H B) \f$
   *
   *  It is defined as
   *  \f$ C_(k+1)  = (I - rho_k ph_k q_k^T) C_k  (I - rho_k q_k p_k^T) + rho_k ph_k p_k^T\f$
   *
   *  For details, we refer to S. Gratton, A. Sartenaer and J. Tshimanga,
   *  SIAM J. Optim.,21(3),912–935,2011 and S. Gurol, PhD Manuscript, 2013.
   *
   *  The solvers represent matrices as objects that implement a "multiply"
   *  method. This class defines objects that apply the preconditioner matrix \f$ P \f$.
   */

// -----------------------------------------------------------------------------

template<typename VECTOR, typename BMATRIX, typename CMATRIX> class QNewtonLMP {
 public:
  explicit QNewtonLMP(const eckit::Configuration &);
  ~QNewtonLMP() {}

  void push(const VECTOR &, const VECTOR &, const VECTOR &, const double &);
  /// Set ObsBias part of the preconditioner to \p Cmat.
  void updateObsBias(std::unique_ptr<CMATRIX> Cmat);
  void update(const BMATRIX & B);

  void multiply(const VECTOR &, VECTOR &) const;

  // Matrix-vector products with the transpose of the preconditioner
  void tmultiply(const VECTOR &, VECTOR &) const;

 private:
  unsigned maxpairs_;
  unsigned maxnewpairs_;
  bool useoldpairs_;
  int maxouter_;
  int update_;

  std::vector<VECTOR> P_;
  std::vector<VECTOR> Ph_;
  std::vector<VECTOR> AP_;
  std::vector<VECTOR> BAP_;
  std::vector<double> rhos_;
  std::vector<unsigned> usedpairIndx_;
  std::unique_ptr<CMATRIX> Cmatrix_;

  std::vector<VECTOR> savedP_;
  std::vector<VECTOR> savedPh_;
  std::vector<VECTOR> savedAP_;
  std::vector<double> savedrhos_;
};

// =============================================================================

template<typename VECTOR, typename BMATRIX, typename CMATRIX>
QNewtonLMP<VECTOR, BMATRIX, CMATRIX>::QNewtonLMP(const eckit::Configuration & conf)
  : maxpairs_(0), maxnewpairs_(0), useoldpairs_(false), maxouter_(0), update_(1)
{
  maxouter_ = conf.getInt("nouter");
  Log::info() << "QNewtonLMP: maxouter : " << maxouter_ << std::endl;

  if (conf.has("preconditioner")) {
    const eckit::LocalConfiguration precond(conf, "preconditioner");
    if (precond.has("maxpairs")) maxpairs_ = precond.getInt("maxpairs");
    if (maxpairs_ > 0) {
      std::string old;
      if (precond.has("useoldpairs")) old = precond.getString("useoldpairs");
      useoldpairs_ = (old == "on" || old == "true");
      if (useoldpairs_ && precond.has("maxnewpairs")) {
          maxnewpairs_ = precond.getInt("maxnewpairs");
          maxnewpairs_ = std::min(maxpairs_, maxnewpairs_);
      } else {
       maxnewpairs_ = maxpairs_;
      }
    }
  }
}

// -----------------------------------------------------------------------------

template<typename VECTOR, typename BMATRIX, typename CMATRIX>
void QNewtonLMP<VECTOR, BMATRIX, CMATRIX>::push(const VECTOR & p, const VECTOR & ph,
                                     const VECTOR & ap, const double & rho) {
  ASSERT(savedP_.size() <= maxnewpairs_);
  if (maxnewpairs_ > 0 && update_ < maxouter_) {
    if (savedP_.size() == maxnewpairs_) {
      savedP_.erase(savedP_.begin());
      savedPh_.erase(savedPh_.begin());
      savedAP_.erase(savedAP_.begin());
      savedrhos_.erase(savedrhos_.begin());
    }
    savedP_.push_back(p);
    savedPh_.push_back(ph);
    savedAP_.push_back(ap);
    savedrhos_.push_back(1.0/rho);
  }
}

// -----------------------------------------------------------------------------

template<typename VECTOR, typename BMATRIX, typename CMATRIX>
void QNewtonLMP<VECTOR, BMATRIX, CMATRIX>::updateObsBias(std::unique_ptr<CMATRIX> Cmat) {
  // Save the preconditioner
  Cmatrix_ = std::move(Cmat);
}

// -----------------------------------------------------------------------------

template<typename VECTOR, typename BMATRIX, typename CMATRIX>
void QNewtonLMP<VECTOR, BMATRIX, CMATRIX>::update(const BMATRIX & Bmat) {
  const unsigned nvec = savedPh_.size();
  ASSERT(savedP_.size() == nvec);
  ASSERT(savedAP_.size() == nvec);
  ASSERT(savedrhos_.size() == nvec);

  if (!useoldpairs_ || nvec >= maxpairs_) {
    Ph_.clear();
    P_.clear();
    AP_.clear();
    BAP_.clear();
    rhos_.clear();
    for (unsigned kiter = 0; kiter < usedpairIndx_.size(); ++kiter) {
      usedpairIndx_[kiter] = 0;
    }
  }

  if (nvec > 0 && update_ < maxouter_) {
    Log::info() << "QNewtonLMP: update " << update_ << ", max = " << maxouter_-1 << std::endl;
    unsigned newpairs = std::min(nvec, maxpairs_);
    unsigned oldpairs = P_.size();
    unsigned rmpairs = 0;
//  First remove pairs we no longer need.
    if (oldpairs + newpairs > maxpairs_) {
      rmpairs = oldpairs + newpairs - maxpairs_;
      for (unsigned jv = 0; jv < rmpairs; ++jv) {
        Ph_.erase(Ph_.begin());
        P_.erase(P_.begin());
        AP_.erase(AP_.begin());
        BAP_.erase(BAP_.begin());
        rhos_.erase(rhos_.begin());
      }
      oldpairs -= rmpairs;
//    Keep information on how many pairs are used at each minimization loop
      unsigned removed = rmpairs;
      for (unsigned kiter = 0; kiter < usedpairIndx_.size(); ++kiter) {
        unsigned rmiter = std::min(usedpairIndx_[kiter], removed);
        usedpairIndx_[kiter] -= rmiter;
        removed -= rmiter;
      }
      ASSERT(removed == 0);
    }
    ASSERT(P_.size() == oldpairs);
    ASSERT(oldpairs + newpairs <= maxpairs_);

//  Add the new pairs.
    for (unsigned jv = nvec - newpairs; jv < nvec; ++jv) {
      P_.push_back(savedP_[jv]);
      Ph_.push_back(savedPh_[jv]);
      AP_.push_back(savedAP_[jv]);
      rhos_.push_back(savedrhos_[jv]);
    // Save B*ap
      VECTOR ww(savedAP_[jv]);
      Bmat.multiply(savedAP_[jv], ww);
      BAP_.push_back(ww);
    }
    ASSERT(P_.size() == oldpairs + newpairs);

    Log::info() << "Number of inner iterations       : " << nvec << std::endl;
    Log::info() << "Number of maximum pairs          : " << maxpairs_ << std::endl;
    Log::info() << "Number of used total pairs       : " << oldpairs + newpairs << std::endl;
    Log::info() << "Number of new pairs              : " << newpairs << std::endl;
    Log::info() << "Number of removed old pairs      : " << rmpairs << std::endl;
    for (unsigned kiter = 0; kiter < usedpairIndx_.size(); ++kiter) {
      Log::info() << "Number of used pairs from outer loop " << kiter + 1
                  << " : "                                   << usedpairIndx_[kiter];
    }
  }

  ++update_;
  savedP_.clear();
  savedPh_.clear();
  savedAP_.clear();
  savedrhos_.clear();
}

// -----------------------------------------------------------------------------
template<typename VECTOR, typename BMATRIX, typename CMATRIX>
void QNewtonLMP<VECTOR, BMATRIX, CMATRIX>::multiply(const VECTOR & a, VECTOR & b) const {
  if (!Cmatrix_) {
    oops::Log::error() << "The VarBC preconditioner matrix is not defined" << std::endl;
    throw eckit::UserError("The VarBC preconditioner matrix is not defined", Here());
  }
  b = a;
  const unsigned nvec = P_.size();
  std::vector<double> etas;
  if (nvec != 0) {
    etas.clear();
    for (int iiter = nvec-1; iiter >= 0; iiter--) {
      double eta = dot_product(b, P_[iiter]);
      eta *= rhos_[iiter];
      etas.push_back(eta);
      b.axpy(-eta, AP_[iiter]);
    }
  }
// For VarBC the obsbias section of the increment is multiplied by the preconditioner
  Cmatrix_->multiply(b, b);
  if (nvec != 0) {
  // Scale b for improved preconditiong of state section of increment.
    b *= dot_product(AP_[nvec-1], AP_[nvec-1])/dot_product(AP_[nvec-1], Ph_[nvec-1]);
    for (unsigned iiter = 0; iiter < nvec; ++iiter) {
      double sigma = dot_product(b, BAP_[iiter]);
      sigma *= rhos_[iiter];
      sigma -= etas[iiter];
      b.axpy(-sigma, Ph_[iiter]);
    }
  }
}

// -----------------------------------------------------------------------------
template<typename VECTOR, typename BMATRIX, typename CMATRIX>
void QNewtonLMP<VECTOR, BMATRIX, CMATRIX>::tmultiply(const VECTOR & a, VECTOR & b) const {
  if (!Cmatrix_) {
    oops::Log::error() << "The VarBC preconditioner matrix is not defined" << std::endl;
    throw eckit::UserError("The VarBC preconditioner matrix is not defined", Here());
  }
  b = a;
  const unsigned nvec = P_.size();
  std::vector<double> etas;
  if (nvec != 0) {
    etas.clear();
    for (int iiter = nvec-1; iiter >= 0; iiter--) {
      double eta = dot_product(b, Ph_[iiter]);
      eta *= rhos_[iiter];
      etas.push_back(eta);
      b.axpy(-eta, BAP_[iiter]);
    }
  }
// For VarBC the obsbias section of the increment is multiplied by the transpose
// of the preconditioner. Note because the preconditioner is curently diagonal
// P^t*b = P*b so we just multiply by the preconditioner.
  Cmatrix_->multiply(b, b);
  if (nvec != 0) {
  // Scale b for improved preconditiong of state section of increment.
    b *= dot_product(AP_[nvec-1], AP_[nvec-1])/dot_product(AP_[nvec-1], Ph_[nvec-1]);
    for (unsigned iiter = 0; iiter < nvec; ++iiter) {
      double sigma = dot_product(b, AP_[iiter]);
      sigma *= rhos_[iiter];
      sigma -= etas[iiter];
      b.axpy(-sigma, P_[iiter]);
    }
  }
}
}  // namespace oops

#endif  // OOPS_ASSIMILATION_QNEWTONLMP_H_
