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
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
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

/// 4D Jb Cost Function
/*!
 * CostJb4D encapsulates the generalized four dimensional Jb term of
 * the 4D-Ens-Var cost function.
 */

template<typename MODEL, typename OBS> class CostJb4D : public CostJbState<MODEL, OBS> {
  typedef ControlIncrement<MODEL, OBS>  CtrlInc_;
  typedef ControlVariable<MODEL, OBS>   CtrlVar_;
  typedef Geometry<MODEL>               Geometry_;
  typedef Increment4D<MODEL>            Increment_;
  typedef State4D<MODEL>                State_;
  typedef JqTerm<MODEL>                 JqTerm_;

 public:
/// Construct \f$ J_b\f$.
  CostJb4D(const std::vector<util::DateTime> &,
           const eckit::Configuration &, const eckit::mpi::Comm &,
           const Geometry_ &, const Variables &);

/// Destructor
  virtual ~CostJb4D() {}

/// Get increment from state (usually first guess).
  void computeIncrement(const CtrlVar_ &, const CtrlVar_ &, const std::shared_ptr<JqTerm_>,
                        CtrlInc_ &) const override;

/// Linearize before the linear computations.
  void linearize(const CtrlVar_ &, const CtrlVar_ &, const Geometry_ &) override;

/// Add Jb gradient.
  void addGradient(const CtrlInc_ &, CtrlInc_ &, CtrlInc_ &) const override;

/// Empty Jq observer.
  JqTermTLAD<MODEL> * initializeJqTLAD() const override {return 0;}

/// Empty TL Jq observer.
  JqTermTLAD<MODEL> * initializeJqTL() const override {return 0;}

/// Empty AD Jq observer.
  JqTermTLAD<MODEL> * initializeJqAD(const CtrlInc_ &) const override {return 0;}

/// Multiply by \f$ B\f$ and \f$ B^{-1}\f$.
  void Bmult(const CtrlInc_ &, CtrlInc_ &) const override;
  void Bminv(const CtrlInc_ &, CtrlInc_ &) const override;

/// Randomize
  void randomize(CtrlInc_ &) const override;

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

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
CostJb4D<MODEL, OBS>::CostJb4D(const std::vector<util::DateTime> & times,
                               const eckit::Configuration & config, const eckit::mpi::Comm & comm,
                               const Geometry_ & geom, const Variables & ctlvars)
  : B_(), bg_(), ctlvars_(ctlvars), resol_(), times_(times), conf_(config, "background error"),
    commTime_(mpi::clone(comm))
{
  bg_.reset(new State_(geom, eckit::LocalConfiguration(config, "background"), commTime_));
  ASSERT(bg_->is_4d());
  ASSERT(bg_->times() == times);
  Log::trace() << "CostJb4D contructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJb4D<MODEL, OBS>::linearize(const CtrlVar_ & xb, const CtrlVar_ & fg,
                                     const Geometry_ & lowres) {
  resol_ = &lowres;
  B_.reset(CovarianceFactory<MODEL>::create(lowres, ctlvars_, conf_, xb.states(), fg.states()));
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJb4D<MODEL, OBS>::computeIncrement(const CtrlVar_ & xb, const CtrlVar_ & fg,
                                       const std::shared_ptr<JqTerm_>, CtrlInc_ & dx) const {
  dx.states().diff(fg.states(), xb.states());
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJb4D<MODEL, OBS>::addGradient(const CtrlInc_ &, CtrlInc_ & grad, CtrlInc_ & gradJb) const {
  grad.states() += gradJb.states();
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJb4D<MODEL, OBS>::Bmult(const CtrlInc_ & dxin, CtrlInc_ & dxout) const {
  B_->multiply(dxin.states(), dxout.states());
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJb4D<MODEL, OBS>::Bminv(const CtrlInc_ & dxin, CtrlInc_ & dxout) const {
  B_->inverseMultiply(dxin.states(), dxout.states());
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJb4D<MODEL, OBS>::randomize(CtrlInc_ & dx) const {
  B_->randomize(dx.states());
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJB4D_H_
