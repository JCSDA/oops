/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_COSTJBJQ_H_
#define OOPS_ASSIMILATION_COSTJBJQ_H_

#include <memory>
#include <vector>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/CostJbState.h"
#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/JqTerm.h"
#include "oops/assimilation/JqTermTLAD.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/util/dot_product.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {
  template<typename MODEL> class ControlIncrement;

// -----------------------------------------------------------------------------

/// Jb + Jq Cost Function
/*!
 * CostJbJq encapsulates the generalized Jb term of the cost weak
 * constraint 4D-Var function (ie Jb+Jq).
 */

template<typename MODEL> class CostJbJq : public CostJbState<MODEL> {
  typedef Increment<MODEL>           Increment_;
  typedef State4D<MODEL>             State4D_;
  typedef Increment4D<MODEL>         Increment4D_;
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef Geometry<MODEL>            Geometry_;

 public:
/// Construct \f$ J_b\f$.
  CostJbJq(const eckit::Configuration &, const Geometry_ &, const Variables &,
           const util::Duration &, const State4D_ &, const bool);

/// Destructor
  virtual ~CostJbJq() {}

/// Finalize \f$ J_q\f$ after the model run.
  JqTerm<MODEL> * initializeJq() const override;

/// Get increment from state (usually first guess).
  void computeIncrement(const State4D_ &, const State4D_ &, Increment4D_ &) const override;

/// Linearize before the linear computations.
  void linearize(const State4D_ &, const Geometry_ &) override;

/// Add Jb gradient.
  void addGradient(const Increment4D_ &, Increment4D_ &, Increment4D_ &) const override;

/// Finalize \f$ J_q\f$ after the TL run.
  JqTermTLAD<MODEL> * initializeJqTL() const override;

/// Initialize \f$ J_q\f$ forcing before the AD run.
  JqTermTLAD<MODEL> * initializeJqAD(const Increment4D_ &) const override;

/// Multiply by \f$ B\f$ and \f$ B^{-1}\f$.
  void Bmult(const Increment4D_ &, Increment4D_ &) const override;
  void Bminv(const Increment4D_ &, Increment4D_ &) const override;

/// Randomize
  void randomize(Increment4D_ &) const override;

/// Create new increment (set to 0).
  unsigned int nstates() const override {return B_.size();}
  Increment_ * newStateIncrement(const unsigned int) const override;

 private:
  boost::ptr_vector< ModelSpaceCovarianceBase<MODEL> > B_;

  const util::Duration subwin_;
  const bool forcing_;
  const Variables controlvars_;
  boost::scoped_ptr<const Geometry_> resol_;
  std::vector<util::DateTime> times_;
};

// =============================================================================

//  Generalized Jb Term of Cost Function
// -----------------------------------------------------------------------------

template<typename MODEL>
CostJbJq<MODEL>::CostJbJq(const eckit::Configuration & config, const Geometry_ & resolouter,
                          const Variables & ctlvars, const util::Duration & len,
                          const State4D_ & xb, const bool forcing)
  : B_(0), subwin_(len), forcing_(forcing), controlvars_(ctlvars), resol_(), times_()
{
  std::vector<eckit::LocalConfiguration> confs;
  config.get("Covariance", confs);
  for (size_t jsub = 0; jsub < confs.size(); ++jsub) {
    B_.push_back(CovarianceFactory<MODEL>::create(confs[jsub], resolouter, ctlvars, xb[jsub]));
  }
  Log::trace() << "CostJbJq contructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbJq<MODEL>::linearize(const State4D_ & fg, const Geometry_ & resolinner) {
  ASSERT(fg.checkStatesNumber(B_.size()));
  resol_.reset(new Geometry_(resolinner));
  times_.clear();
  for (unsigned jsub = 0; jsub < B_.size(); ++jsub) {
    times_.push_back(fg[jsub].validTime());
    B_[jsub].linearize(fg[jsub], *resol_);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbJq<MODEL>::computeIncrement(const State4D_ & xb, const State4D_ & fg,
                                       Increment4D_ & dx) const {
// Compute x_0 - x_b for Jb
  dx[0].diff(fg[0], xb[0]);
  Log::info() << "CostJbJq: x_0 - x_b" << dx[0] << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbJq<MODEL>::addGradient(const Increment4D_ & dxFG, Increment4D_ & grad,
                                  Increment4D_ & gradJb) const {
// Jb from pre-computed gradient
  grad[0] += gradJb[0];

// Compute and add Jq gradient Qi^{-1} ( x_i - M(x_{i-1}) )
  Increment4D_ gg(grad, false);
  for (unsigned jsub = 1; jsub < B_.size(); ++jsub) {
    B_[jsub].inverseMultiply(dxFG[jsub], gg[jsub]);
    grad[jsub] += gg[jsub];
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
JqTerm<MODEL> * CostJbJq<MODEL>::initializeJq() const {
  return new JqTerm<MODEL>(B_.size());
}

// -----------------------------------------------------------------------------

template<typename MODEL>
JqTermTLAD<MODEL> * CostJbJq<MODEL>::initializeJqTL() const {
  JqTermTLAD<MODEL> * jqtl = 0;
  if (!forcing_) jqtl = new JqTermTLAD<MODEL>(B_.size());
  return jqtl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
JqTermTLAD<MODEL> * CostJbJq<MODEL>::initializeJqAD(const Increment4D_ & dx) const {
  JqTermTLAD<MODEL> * jqad = 0;
  if (!forcing_) jqad = new JqTermTLAD<MODEL>(B_.size(), dx);
  return jqad;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbJq<MODEL>::Bmult(const Increment4D_ & dxin, Increment4D_ & dxout) const {
  for (unsigned jsub = 0; jsub < B_.size(); ++jsub) {
    B_[jsub].multiply(dxin[jsub], dxout[jsub]);

    if (jsub == 0) {
      Log::debug() << "CostJbJq:multiply Jb is "
                   << 0.5 * dot_product(dxin[0], dxout[0]) << std::endl;
    } else {
      Log::debug() << "CostJbJq:multiply Jq(" << jsub << ") is "
                   << 0.5 * dot_product(dxin[jsub], dxout[jsub]);
    }
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbJq<MODEL>::randomize(Increment4D_ & dx) const {
  for (unsigned jsub = 0; jsub < B_.size(); ++jsub) {
    B_[jsub].randomize(dx[jsub]);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbJq<MODEL>::Bminv(const Increment4D_ & dxin, Increment4D_ & dxout) const {
  Log::warning() << "*** B inverse might not always exist ***" << std::endl;
  for (unsigned jsub = 0; jsub < B_.size(); ++jsub) {
    B_[jsub].inverseMultiply(dxin[jsub], dxout[jsub]);

    if (jsub == 0) {
      Log::debug() << "CostJbJq:inverseMultiply Jb is "
                   << 0.5 * dot_product(dxin[0], dxout[0]) << std::endl;
    } else {
      Log::debug() << "CostJbJq:inverseMultiply Jq(" << jsub << ") is "
                   << 0.5 * dot_product(dxin[jsub], dxout[jsub]);
    }
  }
  Log::warning() << "*** B inverse might not always exist ***" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> *
CostJbJq<MODEL>::newStateIncrement(const unsigned int isub) const {
  Increment_ * incr = new Increment_(*resol_, controlvars_, times_[isub]);
  return incr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJBJQ_H_
