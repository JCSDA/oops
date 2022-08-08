/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_SPECTRALLMP_H_
#define OOPS_ASSIMILATION_SPECTRALLMP_H_

#include <algorithm>
#include <cmath>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/TriDiagSpectrum.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"

namespace oops {

  /*  The SPECTRAL preconditioner matrix \f$ P \f$ for the system matrix AB = (I + Ht Rinv H B)
   *  It is defined as \f$ P_(k+1)  = P_k  +  U (D^{-1} - I_l) X' \f$
   *  where
   *  \f$ D \f$ is a l-by-l diagonal matrix whose diagonal entries are \f$ lambda_i \f$,
   *  \f$ X = [x_1, x_2, ..., x_l] \f$  and \f$ U = [u_1, u_2, ..., u_l] \f$
   *  are column matrices.
   *
   *  The pairs \f$ (lambda_i,x_i) \f$ define the Ritz pairs of the matrix
   *  \f$ PA \f$ with respect to the subspace \f$ K(PA,Pr0) \f$  and \f$ U = Binv X \f$.
   *   Note that \f$ r0 \f$ is the initial residual.
   *
   *  Note that if RitzPrecond_ is active, it applies the RITZ LMP. For details,
   *  we refer to S. Gratton, A. Sartenaer and J. Tshimanga, SIAM J. Optim.,21(3),912–935,2011
   *  and S. Gurol, PhD Manuscript, 2013.
   */
  /*!
   *  The solvers represent matrices as objects that implement a "multiply"
   *  method. This class defines objects that apply the preconditioner matrix \f$ P \f$.
   */

// -----------------------------------------------------------------------------

template<typename VECTOR, typename CMATRIX> class SpectralLMP {
 public:
  explicit SpectralLMP(const eckit::Configuration &);
  ~SpectralLMP() {}

  /// Set ObsBias part of the preconditioner to \p Cmat.
  void updateObsBias(std::unique_ptr<CMATRIX> Cmat);

  void update(std::vector<std::unique_ptr<VECTOR>> &, std::vector<std::unique_ptr<VECTOR>> &,
              std::vector<std::unique_ptr<VECTOR>> &, std::vector<double> &, std::vector<double> &);

  void multiply(const VECTOR &, VECTOR &) const;

 private:
  unsigned maxpairs_;
  bool useoldpairs_;
  bool RitzPrecond_;
  int maxouter_;
  int update_;

  std::vector<std::unique_ptr<VECTOR>> X_;
  std::vector<std::unique_ptr<VECTOR>> U_;
  std::vector<double> eigvals_;
  std::vector<double> omega_;

  // For RitzPrecond
  std::vector<std::unique_ptr<VECTOR>> Y_;
  std::vector<std::unique_ptr<VECTOR>> S_;
  std::vector<VECTOR> Zlast_;
  std::vector<VECTOR> Zhlast_;
  std::vector<unsigned> usedpairIndx_;
  std::vector<unsigned> zcount;

  // For VarBC preconditioning
  std::unique_ptr<CMATRIX> Cmatrix_;
};

// =============================================================================

template<typename VECTOR, typename CMATRIX>
SpectralLMP<VECTOR, CMATRIX>::SpectralLMP(const eckit::Configuration & conf)
  : maxpairs_(0), useoldpairs_(false), RitzPrecond_(false), maxouter_(0), update_(0)
{
  maxouter_ = conf.getInt("nouter");
  Log::info() << "SpectralLMP: maxouter = " << maxouter_ << std::endl;

  if (conf.has("preconditioner")) {
    const eckit::LocalConfiguration precond(conf, "preconditioner");
    if (precond.has("maxpairs")) maxpairs_ = precond.getInt("maxpairs");

    if (maxpairs_ > 0) {
      Log::info() << "SpectralLMP: maxpairs_ = " << maxpairs_ << std::endl;

      std::string ritz;
      if (precond.has("ritz")) ritz = precond.getString("ritz");
      RitzPrecond_ = (ritz == "on" || ritz == "true");
      Log::info() << "SpectralLMP: RitzPrecond_ = " << RitzPrecond_ << std::endl;

      std::string old;
      if (precond.has("useoldpairs")) old = precond.getString("useoldpairs");
      useoldpairs_ = (old == "on" || old == "true");
      Log::info() << "SpectralLMP: useoldpairs_ = " << useoldpairs_ << std::endl;
    }
  }
}

// -----------------------------------------------------------------------------

template<typename VECTOR, typename CMATRIX>
void SpectralLMP<VECTOR, CMATRIX>::updateObsBias(std::unique_ptr<CMATRIX> Cmat) {
  // Save the preconditioner
  Cmatrix_ = std::move(Cmat);
}

// -----------------------------------------------------------------------------

template<typename VECTOR, typename CMATRIX>
void SpectralLMP<VECTOR, CMATRIX>::update(std::vector<std::unique_ptr<VECTOR>> & Zv,
                                 std::vector<std::unique_ptr<VECTOR>> & Zhl,
                                 std::vector<std::unique_ptr<VECTOR>> & Zl,
                                 std::vector<double> & alphas,
                                 std::vector<double> & betas) {
//  If useoldpairs = false, use only current information
  if (!useoldpairs_) {
    eigvals_.clear();
    omega_.clear();
    X_.clear();
    U_.clear();
    usedpairIndx_.clear();
    Zlast_.clear();
    Zhlast_.clear();
  }

  ++update_;
  const unsigned nvec = Zl.size();
  ASSERT(Zhl.size() == nvec);

  if (nvec > 0) {
    Log::info() << std::endl;
    Log::info() << "SpectralLMP: update " << update_ << ", max = " << maxouter_-1 << std::endl;
    ASSERT(alphas.size() ==  nvec-1);
    ASSERT(betas.size() == nvec-1);
    std::vector<double> evals;
    std::vector< std::vector<double> > evecs;

//  Compute spectrum of tri-diagonal matrix
    TriDiagSpectrum(alphas, betas, evals, evecs);

//  Determine the converged eigenvalues
    std::vector<unsigned> convIndx;
    unsigned oldpairIndx = 0;
    for (unsigned jiter = 0; jiter < nvec-1; ++jiter) {
      double erritz = evecs[jiter][nvec-2]*betas[nvec-2];
      double lambda = evals[jiter];
//    Add new converged eigenvalues
      if (std::abs(erritz) <= 0.001*lambda) {
        convIndx.push_back(jiter);
        eigvals_.push_back(lambda);
        omega_.push_back(erritz/lambda);
//      If maximum number of pairs are exceeded, remove old information from
//      the beginning
        if (convIndx.size() > maxpairs_) {
          convIndx.erase(convIndx.begin());
        }
        if (eigvals_.size() > maxpairs_) {
          eigvals_.erase(eigvals_.begin());
          omega_.erase(omega_.begin());
          oldpairIndx += 1;
        }
      }
    }
//  Keep information on how many pairs are used at each minimization loop
    usedpairIndx_.push_back(convIndx.size());
    unsigned reduc = oldpairIndx;
    for (unsigned kiter = 0; kiter < usedpairIndx_.size()-1; ++kiter) {
      if (reduc <= usedpairIndx_[kiter]) {
        usedpairIndx_[kiter] -= reduc;
        break;
      } else {
        reduc -= usedpairIndx_[kiter];
        usedpairIndx_[kiter] = 0;
      }
    }
    unsigned counter = 0;
    for (unsigned kiter = 0; kiter < usedpairIndx_.size(); ++kiter) {
      counter += usedpairIndx_[kiter];
    }
    ASSERT(counter <=  maxpairs_);
    Log::info() << "SpectralLMP: Number of inner iterations    : " << Zl.size()-1 << std::endl;
    Log::info() << "SpectralLMP: Number of maximum pairs       : " << maxpairs_ << std::endl;
    Log::info() << "SpectralLMP: Number of used total pairs    : " << counter << std::endl;
    Log::info() << "SpectralLMP: Number of new converged pairs : " << convIndx.size() << std::endl;
    Log::info() << "SpectralLMP: Number of removed old pairs   : " << oldpairIndx << std::endl;
    for (unsigned kiter = 0; kiter < usedpairIndx_.size(); ++kiter) {
      Log::info() << "SpectralLMP: Number of used pairs from outer loop " << kiter + 1
                << " : " << usedpairIndx_[kiter] << std::endl;
    }

//  Calculate/Update the matrix X and U

    if (X_.size() != 0) {
//    Remove the first oldpairIndx elements
      unsigned minsize = std::min(oldpairIndx, maxpairs_);
      unsigned xsize = X_.size();
      X_.erase(X_.begin(), X_.begin() + std::min(xsize, minsize));
      U_.erase(U_.begin(), U_.begin() + std::min(xsize, minsize));
    }
    for (unsigned jiter = 0; jiter < convIndx.size(); ++jiter) {
      std::unique_ptr<VECTOR> ww(new VECTOR(*Zl[0], false));
      for (unsigned iiter = 0; iiter < nvec-1; ++iiter) {
        ww->axpy(evecs[convIndx[jiter]][iiter], *Zl[iiter]);
      }
//    Add new information
      X_.emplace_back(std::move(ww));
    }
    for (unsigned jiter = 0; jiter < convIndx.size(); ++jiter) {
      std::unique_ptr<VECTOR> ww(new VECTOR(*Zl[0], false));
      for (unsigned iiter = 0; iiter < nvec-1; ++iiter) {
        ww->axpy(evecs[convIndx[jiter]][iiter], *Zhl[iiter]);
      }
//    Add new information
      U_.emplace_back(std::move(ww));
    }

    if (RitzPrecond_) {
      Zlast_.push_back(*Zl[nvec-1]);
      Zhlast_.push_back(*Zhl[nvec-1]);

//    Calculate the matrix Y = [U1*omega1, ..., Uk*omegak]
      Y_.clear();
      zcount.clear();
      unsigned kk = 0;
      for (unsigned kiter = 0; kiter < usedpairIndx_.size(); ++kiter) {
        if (usedpairIndx_[kiter] != 0) {
          zcount.push_back(kiter);
          std::unique_ptr<VECTOR> ww(new VECTOR(*Zl[0], false));
          for (unsigned jiter = 0; jiter < usedpairIndx_[kiter]; ++jiter) {
            ww->axpy(omega_[jiter + kk], *U_[jiter + kk]);
          }
          kk += usedpairIndx_[kiter];
          Y_.emplace_back(std::move(ww));
        }
      }

//    Calculate the matrix S = [X1*omega1, ..., Xk*omegak]
      S_.clear();
      kk = 0;
      for (unsigned kiter = 0; kiter < usedpairIndx_.size(); ++kiter) {
        if (usedpairIndx_[kiter] != 0) {
          std::unique_ptr<VECTOR> ww(new VECTOR(*Zl[0], false));
          for (unsigned jiter = 0; jiter < usedpairIndx_[kiter]; ++jiter) {
            ww->axpy(omega_[jiter + kk], *X_[jiter + kk]);
          }
          kk += usedpairIndx_[kiter];
          S_.emplace_back(std::move(ww));
        }
      }
    }

//  Release remaining memory
    Zv.clear();
    Zl.clear();
    Zhl.clear();
    alphas.clear();
    betas.clear();
  }
}

// -----------------------------------------------------------------------------

template<typename VECTOR, typename CMATRIX>
void SpectralLMP<VECTOR, CMATRIX>::multiply(const VECTOR & a, VECTOR & b) const {
  b = a;  // P_0 = I
  // For VarBC the obsbias section of the increment is multiplied by the preconditioner
  if (!Cmatrix_) {
    oops::Log::error() << "The VarBC preconditioner matrix is not defined" << std::endl;
    throw eckit::UserError("The VarBC preconditioner matrix is not defined", Here());
  }
  Cmatrix_->multiply(a, b);

  for (unsigned iiter = 0; iiter < eigvals_.size(); ++iiter) {
    double zeval = std::min(10.0, eigvals_[iiter]);
//    double zeval = eigvals_[iiter];
    double zz = 1.0/zeval - 1.0;
    b.axpy(zz*dot_product(a, *X_[iiter]), *U_[iiter]);
  }

  if (RitzPrecond_ && eigvals_.size() != 0) {
//  Sk'*a
    std::vector<double> sta;
    sta.clear();
    for (unsigned iiter = 0; iiter < S_.size(); ++iiter) {
      sta.push_back(dot_product(*S_[iiter], a));
    }

//  Yk (sta(k) - Zlast' a)
    for (unsigned kiter = 0; kiter < S_.size(); ++kiter) {
      for (unsigned iiter = 0; iiter < eigvals_.size(); ++iiter) {
        double wxap = sta[kiter] - dot_product(a, Zlast_[zcount[kiter]]);
        b.axpy(wxap, *Y_[kiter]);
      }
    }

//  -Zhlast sta
    for (unsigned kiter = 0; kiter < S_.size(); ++kiter) {
      b.axpy(-sta[kiter], Zhlast_[zcount[kiter]]);
    }
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_SPECTRALLMP_H_
