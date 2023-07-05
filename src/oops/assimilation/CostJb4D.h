/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2021-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_COSTJB4D_H_
#define OOPS_ASSIMILATION_COSTJB4D_H_

#include <memory>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/CostJbState.h"
#include "oops/assimilation/GMRESR.h"
#include "oops/base/Geometry.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/Increment4D.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace oops {
  template<typename MODEL> class JqTerm;
  template<typename MODEL> class JqTermTLAD;

// -----------------------------------------------------------------------------

template<typename MODEL> class Bwrapper {
  typedef Increment4D<MODEL>               Increment4D_;
  typedef ModelSpaceCovarianceBase<MODEL>  Covariance_;

 public:
  explicit Bwrapper(const Covariance_ & B) : B_(&B) {}
  ~Bwrapper() {}

  void multiply(const Increment4D_ & dxin, Increment4D_ & dxout) const {
    B_->multiply(dxin[0], dxout[0]);
  }

 private:
  const Covariance_ * B_;
};

// -----------------------------------------------------------------------------

/// 4D Jb Cost Function
/*!
 * CostJb4D encapsulates the generalized four dimensional Jb term of
 * the 4D-Ens-Var cost function.
 */

template<typename MODEL> class CostJb4D : public CostJbState<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment4D<MODEL>         Increment_;
  typedef State4D<MODEL>             State_;
  typedef JqTerm<MODEL>              JqTerm_;

 public:
/// Construct \f$ J_b\f$.
  CostJb4D(const std::vector<util::DateTime> &,
           const eckit::Configuration &, const eckit::mpi::Comm &,
           const Geometry_ &, const Variables &);

/// Destructor
  virtual ~CostJb4D() {}

/// Get increment from state (usually first guess).
  void computeIncrement(const State_ &, const State_ &, const std::shared_ptr<JqTerm_>,
                        Increment_ &) const override;

/// Linearize before the linear computations.
  void linearize(const State_ &, const State_ &, const Geometry_ &) override;

/// Add Jb gradient.
  void addGradient(const Increment_ &, Increment_ &, Increment_ &) const override;

/// Empty Jq observer.
  JqTermTLAD<MODEL> * initializeJqTLAD() const override {return 0;}

/// Empty TL Jq observer.
  JqTermTLAD<MODEL> * initializeJqTL() const override {return 0;}

/// Empty AD Jq observer.
  JqTermTLAD<MODEL> * initializeJqAD(const Increment_ &) const override {return 0;}

/// Multiply by \f$ B\f$ and \f$ B^{-1}\f$.
  void Bmult(const Increment_ &, Increment_ &) const override;
  void Bminv(const Increment_ &, Increment_ &) const override;

/// Randomize
  void randomize(Increment_ &) const override;

/// Accessors to data for constructing a new increment.
  const Geometry_ & geometry() const override {return *resol_;}
  const Variables & variables() const override {return ctlvars_;}
  const std::vector<util::DateTime> & times() const override {return times_;}
  const eckit::mpi::Comm & comm() const override {return commTime_;}
  std::shared_ptr<State_> background() const override {return bg_;}

 private:
  std::unique_ptr<ModelSpaceCovarianceBase<MODEL>> B_;
  std::shared_ptr<State_> bg_;
  const Variables ctlvars_;
  const Geometry_ * resol_;
  std::vector<util::DateTime> times_;
  const eckit::LocalConfiguration conf_;
  const eckit::mpi::Comm & commTime_;
};

// =============================================================================

//  Generalized Jb Term of Cost Function
// -----------------------------------------------------------------------------

template<typename MODEL>
CostJb4D<MODEL>::CostJb4D(const std::vector<util::DateTime> & times,
                          const eckit::Configuration & config, const eckit::mpi::Comm & comm,
                          const Geometry_ & geom, const Variables & ctlvars)
  : B_(), bg_(), ctlvars_(ctlvars), resol_(), times_(times), conf_(config, "background error"),
    commTime_(mpi::clone(comm))
{
  bg_.reset(new State_(times, commTime_)),
  bg_->read(geom, eckit::LocalConfiguration(config, "background"));
  ASSERT(bg_->is_4d());
  Log::trace() << "CostJb4D contructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::linearize(const State_ & xb, const State_ & fg, const Geometry_ & lowres) {
  resol_ = &lowres;
  B_.reset(CovarianceFactory<MODEL>::create(lowres, ctlvars_, conf_, xb[0], fg[0]));
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::computeIncrement(const State_ & xb, const State_ & fg,
                                       const std::shared_ptr<JqTerm_>, Increment_ & dx) const {
  dx.diff(fg, xb);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::addGradient(const Increment_ &, Increment_ & grad,
                                  Increment_ & gradJb) const {
  grad += gradJb;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::Bmult(const Increment_ & dxin, Increment_ & dxout) const {
  B_->multiply(dxin[0], dxout[0]);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::Bminv(const Increment_ & dxin, Increment_ & dxout) const {
//  B_->inverseMultiply(dxin[0], dxout[0]);
  Bwrapper<MODEL> BW(*B_);
  IdentityMatrix<Increment_> Id;
  dxout.zero();
  GMRESR(dxout, dxin, BW, Id, 10, 1.0e-3);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::randomize(Increment_ & dx) const {
  B_->randomize(dx[0]);
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJB4D_H_
