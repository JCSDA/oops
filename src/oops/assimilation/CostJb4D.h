/*
 * (C) Copyright 2009-2016 ECMWF.
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
#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/ModelSpaceCovariance4DBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Logger.h"

namespace oops {
  template<typename MODEL> class JqTerm;
  template<typename MODEL> class JqTermTLAD;

// -----------------------------------------------------------------------------

/// 4D Jb Cost Function
/*!
 * CostJb4D encapsulates the generalized four dimensional Jb term of
 * the 4D-Ens-Var cost function.
 */

template<typename MODEL> class CostJb4D : public CostJbState<MODEL> {
  typedef Increment<MODEL>           Increment_;
  typedef State4D<MODEL>             State4D_;
  typedef Increment4D<MODEL>         Increment4D_;
  typedef Geometry<MODEL>            Geometry_;

 public:
/// Construct \f$ J_b\f$.
  CostJb4D(const eckit::Configuration &, const Geometry_ &, const Variables &, const State4D_ &);

/// Destructor
  virtual ~CostJb4D() {}

/// Empty Jq observer.
  JqTerm<MODEL> * initializeJq() const override {return 0;}
  JqTermTLAD<MODEL> * initializeJqTLAD() const override {return 0;}

/// Get increment from state (usually first guess).
  void computeIncrement(const State4D_ &, const State4D_ &, Increment4D_ &) const override;

/// Linearize before the linear computations.
  void linearize(const State4D_ &, const Geometry_ &) override;

/// Add Jb gradient.
  void addGradient(const Increment4D_ &, Increment4D_ &, Increment4D_ &) const override;

/// Empty TL Jq observer.
  JqTermTLAD<MODEL> * initializeJqTL() const override {return 0;}

/// Empty AD Jq observer.
  JqTermTLAD<MODEL> * initializeJqAD(const Increment4D_ &) const override {return 0;}

/// Multiply by \f$ B\f$ and \f$ B^{-1}\f$.
  void Bmult(const Increment4D_ &, Increment4D_ &) const override;
  void Bminv(const Increment4D_ &, Increment4D_ &) const override;

/// Randomize
  void randomize(Increment4D_ &) const override;

/// Create new increment (set to 0).
  unsigned int nstates() const override {return xb_.size();}
  Increment_ * newStateIncrement(const unsigned int) const override;

 private:
  const State4D_ & xb_;
  std::unique_ptr<ModelSpaceCovariance4DBase<MODEL> > B_;
  const Variables ctlvars_;
  std::unique_ptr<const Geometry_> resol_;
  std::vector<util::DateTime> times_;
  const eckit::LocalConfiguration conf_;
};

// =============================================================================

//  Generalized Jb Term of Cost Function
// -----------------------------------------------------------------------------

template<typename MODEL>
CostJb4D<MODEL>::CostJb4D(const eckit::Configuration & config, const Geometry_ &,
                          const Variables & ctlvars, const State4D_ & xb)
  : xb_(xb), B_(), ctlvars_(ctlvars), resol_(), times_(), conf_(config, "background error")
{
  Log::trace() << "CostJb4D contructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::linearize(const State4D_ & fg, const Geometry_ & lowres) {
  ASSERT(fg.checkStatesNumber(xb_.size()));
  resol_.reset(new Geometry_(lowres));
  times_ = fg.validTimes();
  B_.reset(Covariance4DFactory<MODEL>::create(conf_, lowres, ctlvars_, xb_, fg));
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::computeIncrement(const State4D_ & xb, const State4D_ & fg,
                                       Increment4D_ & dx) const {
  dx.diff(fg, xb);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::addGradient(const Increment4D_ & dxFG, Increment4D_ & grad,
                                  Increment4D_ & gradJb) const {
  grad += gradJb;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::Bmult(const Increment4D_ & dxin, Increment4D_ & dxout) const {
  B_->multiply(dxin, dxout);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::Bminv(const Increment4D_ & dxin, Increment4D_ & dxout) const {
  B_->inverseMultiply(dxin, dxout);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::randomize(Increment4D_ & dx) const {
  B_->randomize(dx);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> *
CostJb4D<MODEL>::newStateIncrement(const unsigned int isub) const {
  Increment_ * incr = new Increment_(*resol_, ctlvars_, times_[isub]);
  return incr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJB4D_H_
