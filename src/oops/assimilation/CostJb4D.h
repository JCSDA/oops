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
  boost::ptr_vector< ModelSpaceCovarianceBase<MODEL> > B_;
  const Variables controlvars_;
  boost::scoped_ptr<const Geometry_> resol_;
  std::vector<util::DateTime> times_;
};

// =============================================================================

//  Generalized Jb Term of Cost Function
// -----------------------------------------------------------------------------

template<typename MODEL>
CostJb4D<MODEL>::CostJb4D(const eckit::Configuration & config, const Geometry_ & resolouter,
                          const Variables & ctlvars, const util::Duration &, const State4D_ & xb)
  : B_(), controlvars_(ctlvars), resol_(), times_()
{
// Create one row of blocks of the whole BMatrix, one object for each
// subwindow. Each object stands for all blocks in the same column of the B Matrix.
// It can be from any concrete class of ModelSpaceCovarianceBase,
// according to the configuration file.
  const eckit::LocalConfiguration covar(config, "Covariance");
  std::vector<eckit::LocalConfiguration> confs;
  covar.get("covariance_time", confs);

  for (size_t jsub = 0; jsub < confs.size(); ++jsub) {
    B_.push_back(CovarianceFactory<MODEL>::create(confs[jsub], resolouter, ctlvars, xb[jsub]));
  }

  Log::trace() << "CostJb4D contructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::linearize(const State4D_ & fg, const Geometry_ & resolinner) {
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
  Increment_ * incr = new Increment_(*resol_, controlvars_, times_[isub]);
  return incr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJB4D_H_
