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
#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/CostJbState.h"
#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {
  template<typename MODEL> class ControlIncrement;
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
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef Geometry<MODEL>            Geometry_;

 public:
/// Construct \f$ J_b\f$.
  CostJb4D(const eckit::Configuration &, const Geometry_ &, const Variables &,
           const util::Duration &, const State4D_ &);

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
  unsigned int nstates() const override {return B_.size();}
  Increment_ * newStateIncrement(const unsigned int) const override;

 private:
  const State4D_ & xb_;
  boost::ptr_vector< ModelSpaceCovarianceBase<MODEL> > B_;
  const Variables ctlvars_;
  boost::scoped_ptr<const Geometry_> resol_;
  std::vector<util::DateTime> times_;
  const eckit::LocalConfiguration conf_;
};

// =============================================================================

//  Generalized Jb Term of Cost Function
// -----------------------------------------------------------------------------

template<typename MODEL>
CostJb4D<MODEL>::CostJb4D(const eckit::Configuration & config, const Geometry_ &,
                          const Variables & ctlvars, const util::Duration &, const State4D_ & xb)
  : xb_(xb), B_(), ctlvars_(ctlvars), resol_(), times_(), conf_(config, "Covariance")
{
  Log::trace() << "CostJb4D contructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::linearize(const State4D_ & fg, const Geometry_ & lowres) {
  ASSERT(fg.checkStatesNumber(xb_.size()));
  resol_.reset(new Geometry_(lowres));
  times_.clear();
  B_.clear();
  std::vector<eckit::LocalConfiguration> confs;
  conf_.get("covariance_time", confs);
  for (unsigned jsub = 0; jsub < fg.size(); ++jsub) {
    B_.push_back(CovarianceFactory<MODEL>::create(confs[jsub], lowres, ctlvars_,
                                                  xb_[jsub], fg[jsub]));
    times_.push_back(fg[jsub].validTime());
  }
  ASSERT(fg.checkStatesNumber(B_.size()));
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
  for (unsigned k1 = 0; k1 < B_.size(); ++k1) {
    // Apply the line k1 of the whole B matrix to the StateIncrement dxin
    // Result is the part of dxout at time k1
    dxout[k1].zero();
    Increment_ dout(dxout[k1]);
    for (unsigned k2 = 0; k2 < B_.size(); ++k2) {
      // Apply to increment at time k2 the block B_k1k2 of the whole B matrix.
      // We need an object B which is also at time k2
      // Result "dout" is an increment at time k1
      B_[k2].multiply(dxin[k2], dout);
      dxout[k1] += dout;
    }
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::Bminv(const Increment4D_ & dxin, Increment4D_ & dxout) const {
  Log::warning() << "*** B inverse might not always exist ***" << std::endl;
  for (unsigned jsub = 0; jsub < B_.size(); ++jsub) {
    B_[jsub].inverseMultiply(dxin[jsub], dxout[jsub]);
  }
  Log::warning() << "*** B inverse might not always exist ***" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::randomize(Increment4D_ & dx) const {
  for (unsigned jsub = 0; jsub < B_.size(); ++jsub) {
    B_[jsub].randomize(dx[jsub]);
  }
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
