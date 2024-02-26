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

#ifndef OOPS_ASSIMILATION_COSTJB3D_H_
#define OOPS_ASSIMILATION_COSTJB3D_H_

#include <memory>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostJbState.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment4D.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/State4D.h"
#include "oops/base/Variables.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {
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

template<typename MODEL, typename OBS> class CostJb3D : public CostJbState<MODEL, OBS> {
  typedef ControlIncrement<MODEL, OBS>  CtrlInc_;
  typedef ControlVariable<MODEL, OBS>   CtrlVar_;
  typedef Geometry<MODEL>               Geometry_;
  typedef State4D<MODEL>                State_;
  typedef JqTerm<MODEL>                 JqTerm_;

 public:
/// Construct \f$ J_b\f$.
  CostJb3D(const util::DateTime &, const eckit::Configuration &,
           const Geometry_ &, const Variables &);

/// Destructor
  virtual ~CostJb3D() {}

/// Get increment from state (usually first guess).
  void computeIncrement(const CtrlVar_ &, const CtrlVar_ &, const std::shared_ptr<JqTerm_>,
                        CtrlInc_ &) const override;

/// Linearize before the linear computations.
  void linearize(const CtrlVar_ &, const CtrlVar_ &, const Geometry_ &) override;

/// Add Jb gradient.
  void addGradient(const CtrlInc_ &, CtrlInc_ &, CtrlInc_ &) const override;

  void setTime(const CtrlVar_ & xx) override {time_[0] = xx.state(0).validTime();}

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
  const std::vector<util::DateTime> & times() const override {return time_;}
  const eckit::mpi::Comm & comm() const override {return oops::mpi::myself();}
  std::shared_ptr<State_> background() const override {return bg_;}

 private:
  std::unique_ptr<ModelSpaceCovarianceBase<MODEL>> B_;
  std::shared_ptr<State_> bg_;
  const Variables ctlvars_;
  const Geometry_ * resol_;
  std::vector<util::DateTime> time_;
  const eckit::LocalConfiguration conf_;
};

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
CostJb3D<MODEL, OBS>::CostJb3D(const util::DateTime & time, const eckit::Configuration & config,
                          const Geometry_ & geom, const Variables & ctlvars)
  : B_(), bg_(), ctlvars_(ctlvars), resol_(), time_(1), conf_(config, "background error")
{
  bg_.reset(new State_(geom, eckit::LocalConfiguration(config, "background"), oops::mpi::myself()));
  ASSERT(bg_->is_3d());
  ASSERT(bg_->validTimes()[0] == time);
  Log::trace() << "CostJb3D constructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJb3D<MODEL, OBS>::linearize(const CtrlVar_ & xb, const CtrlVar_ & fg,
                                     const Geometry_ & lowres) {
  Log::trace() << "CostJb3D:linearize start" << std::endl;
  resol_ = &lowres;
  time_[0] = xb.state(0).validTime();  // not earlier because of FGAT
  B_.reset(CovarianceFactory<MODEL>::create(lowres, ctlvars_, conf_, xb.states(), fg.states()));
  Log::trace() << "CostJb3D:linearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJb3D<MODEL, OBS>::computeIncrement(const CtrlVar_ & xb, const CtrlVar_ & fg,
                                       const std::shared_ptr<JqTerm_>, CtrlInc_ & dx) const {
  dx.states().diff(fg.states(), xb.states());
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJb3D<MODEL, OBS>::addGradient(const CtrlInc_ & dxFG, CtrlInc_ & grad,
                                  CtrlInc_ & gradJb) const {
  grad.states() += gradJb.states();
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJb3D<MODEL, OBS>::Bmult(const CtrlInc_ & dxin, CtrlInc_ & dxout) const {
  B_->multiply(dxin.states(), dxout.states());
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJb3D<MODEL, OBS>::Bminv(const CtrlInc_ & dxin, CtrlInc_ & dxout) const {
  B_->inverseMultiply(dxin.states(), dxout.states());
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJb3D<MODEL, OBS>::randomize(CtrlInc_ & dx) const {
  B_->randomize(dx.states());
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJB3D_H_
