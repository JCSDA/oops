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

template<typename MODEL> class CostJb3D : public CostJbState<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment4D<MODEL>         Increment_;
  typedef State4D<MODEL>             State_;
  typedef JqTerm<MODEL>              JqTerm_;

 public:
/// Construct \f$ J_b\f$.
  CostJb3D(const util::DateTime &, const eckit::Configuration &,
           const Geometry_ &, const Variables &);

/// Destructor
  virtual ~CostJb3D() {}

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
  const Geometry_ & geometry() const override;
  const Variables & variables() const override;
  const std::vector<util::DateTime> & times() const override;
  const eckit::mpi::Comm & comm() const override {return oops::mpi::myself();}
  std::shared_ptr<State_> background() const override {return bg_;}

 private:
  std::unique_ptr<ModelSpaceCovarianceBase<MODEL>> B_;
  std::shared_ptr<State_> bg_;
  const Variables controlvars_;
  const Geometry_ * resol_;
  std::vector<util::DateTime> time_;
  const eckit::LocalConfiguration conf_;
};

// =============================================================================

//  Jb Term of Cost Function
// -----------------------------------------------------------------------------

template<typename MODEL>
CostJb3D<MODEL>::CostJb3D(const util::DateTime & time, const eckit::Configuration & config,
                          const Geometry_ & geom, const Variables & ctlvars)
  : B_(), bg_(), controlvars_(ctlvars), resol_(), time_(1), conf_(config, "background error")
{
  std::vector<util::DateTime> times({time});
  bg_.reset(new State_(times, oops::mpi::myself()));
  bg_->read(geom, eckit::LocalConfiguration(config, "background"));
  ASSERT(bg_->is_3d());
  Log::trace() << "CostJb3D constructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb3D<MODEL>::linearize(const State_ & xb, const State_ & fg, const Geometry_ & lowres) {
  Log::trace() << "CostJb3D:linearize start" << std::endl;
  resol_ = &lowres;
  time_[0] = xb[0].validTime();  // not earlier because of FGAT
  B_.reset(CovarianceFactory<MODEL>::create(lowres, controlvars_, conf_, xb[0], fg[0]));
  Log::trace() << "CostJb3D:linearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb3D<MODEL>::computeIncrement(const State_ & xb, const State_ & fg,
                                       const std::shared_ptr<JqTerm_>, Increment_ & dx) const {
  dx.diff(fg, xb);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb3D<MODEL>::addGradient(const Increment_ & dxFG, Increment_ & grad,
                                  Increment_ & gradJb) const {
  grad += gradJb;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb3D<MODEL>::Bmult(const Increment_ & dxin, Increment_ & dxout) const {
  B_->multiply(dxin[0], dxout[0]);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb3D<MODEL>::Bminv(const Increment_ & dxin, Increment_ & dxout) const {
  B_->inverseMultiply(dxin[0], dxout[0]);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJb3D<MODEL>::randomize(Increment_ & dx) const {
  B_->randomize(dx[0]);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
const Geometry<MODEL> & CostJb3D<MODEL>::geometry() const {
  return *resol_;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
const Variables & CostJb3D<MODEL>::variables() const {
  return controlvars_;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
const std::vector<util::DateTime> & CostJb3D<MODEL>::times() const {
  return time_;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJB3D_H_
