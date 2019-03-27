/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_COSTJB3D_H_
#define OOPS_ASSIMILATION_COSTJB3D_H_

#include <memory>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/CostJbState.h"
#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {
  template<typename MODEL> class ControlIncrement;
  template<typename MODEL> class JqTerm;
  template<typename MODEL> class JqTermTLAD;

// -----------------------------------------------------------------------------

/// Jb Cost Function
/*!
 * The CostJb3D encapsulates the Jb term of the cost function for a
 * 3 dimensional background.
 *
 * This class is not really necessary since it is only a special
 * case of the more general CostJbJq weak constraint term
 * with one sub-window. It is provided for readability.
 */

template<typename MODEL> class CostJb3D : public CostJbState<MODEL> {
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;
  typedef State4D<MODEL>             State4D_;
  typedef Increment4D<MODEL>         Increment4D_;
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef Geometry<MODEL>            Geometry_;

 public:
/// Construct \f$ J_b\f$.
  CostJb3D(const eckit::Configuration &, const Geometry_ &, const Variables &,
           const util::Duration &, const State_ &);

/// Destructor
  virtual ~CostJb3D() {}

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
  unsigned int nstates() const override {return 1;}
  Increment_ * newStateIncrement(const unsigned int) const override;

 private:
  const State_ & xb_;
  boost::scoped_ptr< ModelSpaceCovarianceBase<MODEL> > B_;
  const util::Duration winLength_;
  const Variables controlvars_;
  boost::scoped_ptr<const Geometry_> resol_;
  boost::scoped_ptr<const util::DateTime> time_;
  const eckit::LocalConfiguration conf_;
};

// =============================================================================

//  Jb Term of Cost Function
// -----------------------------------------------------------------------------

template<typename MODEL>
CostJb3D<MODEL>::CostJb3D(const eckit::Configuration & config, const Geometry_ &,
                          const Variables & ctlvars, const util::Duration & len,
                          const State_ & xb)
  : xb_(xb), B_(), winLength_(len), controlvars_(ctlvars), resol_(), time_(),
    conf_(config, "Covariance")
{
  Log::trace() << "CostJb3D constructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb3D<MODEL>::linearize(const State4D_ & fg, const Geometry_ & lowres) {
  ASSERT(fg.checkStatesNumber(1));
  resol_.reset(new Geometry_(lowres));
  time_.reset(new util::DateTime(fg[0].validTime()));
  B_.reset(CovarianceFactory<MODEL>::create(conf_, lowres, controlvars_, xb_, fg[0]));
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb3D<MODEL>::computeIncrement(const State4D_ & xb, const State4D_ & fg,
                                       Increment4D_ & dx) const {
  dx.diff(fg, xb);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb3D<MODEL>::addGradient(const Increment4D_ & dxFG, Increment4D_ & grad,
                                  Increment4D_ & gradJb) const {
  grad += gradJb;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb3D<MODEL>::Bmult(const Increment4D_ & dxin, Increment4D_ & dxout) const {
  B_->multiply(dxin[0], dxout[0]);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb3D<MODEL>::Bminv(const Increment4D_ & dxin, Increment4D_ & dxout) const {
  B_->inverseMultiply(dxin[0], dxout[0]);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb3D<MODEL>::randomize(Increment4D_ & dx) const {
  B_->randomize(dx[0]);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> *
CostJb3D<MODEL>::newStateIncrement(const unsigned int) const {
  Increment_ * incr = new Increment_(*resol_, controlvars_, *time_);
  return incr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJB3D_H_
