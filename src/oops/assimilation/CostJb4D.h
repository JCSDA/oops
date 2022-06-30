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
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace oops {
  template<typename MODEL> class JqTermTLAD;

// -----------------------------------------------------------------------------

/// 4D Jb Cost Function
/*!
 * CostJb4D encapsulates the generalized four dimensional Jb term of
 * the 4D-Ens-Var cost function.
 */

template<typename MODEL> class CostJb4D : public CostJbState<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
/// Construct \f$ J_b\f$.
  CostJb4D(const eckit::Configuration &, const eckit::mpi::Comm &,
           const Geometry_ &, const Variables &, const State_ &);

/// Destructor
  virtual ~CostJb4D() {}

/// Get increment from state (usually first guess).
  void computeIncrement(const State_ &, const State_ &, const State_ &,
                        Increment_ &) const override;

/// Linearize before the linear computations.
  void linearize(const State_ &, const Geometry_ &) override;

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

/// Create new increment (set to 0).
  Increment_ * newStateIncrement() const override;

 private:
  const State_ & xb_;
  std::unique_ptr<ModelSpaceCovarianceBase<MODEL> > B_;
  const Variables ctlvars_;
  const Geometry_ * resol_;
  util::DateTime time_;
  const eckit::LocalConfiguration conf_;
  const eckit::mpi::Comm & commTime_;
};

// =============================================================================

//  Generalized Jb Term of Cost Function
// -----------------------------------------------------------------------------

template<typename MODEL>
CostJb4D<MODEL>::CostJb4D(const eckit::Configuration & config, const eckit::mpi::Comm & comm,
                          const Geometry_ &, const Variables & ctlvars, const State_ & xb)
  : xb_(xb), B_(), ctlvars_(ctlvars), resol_(), time_(xb.validTime()),
    conf_(config, "background error"), commTime_(comm)
{
  Log::trace() << "CostJb4D contructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::linearize(const State_ & fg, const Geometry_ & lowres) {
  resol_ = &lowres;
  B_.reset(CovarianceFactory<MODEL>::create(lowres, ctlvars_, conf_, xb_, fg));
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::computeIncrement(const State_ & xb, const State_ & fg, const State_ &,
                                       Increment_ & dx) const {
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
  B_->multiply(dxin, dxout);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::Bminv(const Increment_ & dxin, Increment_ & dxout) const {
  B_->inverseMultiply(dxin, dxout);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb4D<MODEL>::randomize(Increment_ & dx) const {
  B_->randomize(dx);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> *
CostJb4D<MODEL>::newStateIncrement() const {
  Increment_ * incr = new Increment_(*resol_, ctlvars_, time_);
  return incr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJB4D_H_
