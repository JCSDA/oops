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

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/CostJbState.h"
#include "oops/assimilation/JqTermTLAD.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Jb + Jq Cost Function
/*!
 * CostJbJq encapsulates the generalized Jb term of the cost weak
 * constraint 4D-Var function (ie Jb+Jq).
 */

template<typename MODEL> class CostJbJq : public CostJbState<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
/// Construct \f$ J_b\f$.
  CostJbJq(const eckit::Configuration &, const eckit::mpi::Comm &,
           const Geometry_ &, const Variables &, const State_ &);

/// Destructor
  virtual ~CostJbJq() {}

/// Get increment from state (usually first guess).
  void computeIncrement(const State_ &, const State_ &, const State_ &,
                        Increment_ &) const override;

/// Linearize before the linear computations.
  void linearize(const State_ &, const Geometry_ &) override;

/// Add Jb gradient.
  void addGradient(const Increment_ &, Increment_ &, Increment_ &) const override;

/// Finalize \f$ J_q\f$ after the model run.
  JqTermTLAD<MODEL> * initializeJqTLAD() const override;

/// Finalize \f$ J_q\f$ after the TL run.
  JqTermTLAD<MODEL> * initializeJqTL() const override;

/// Initialize \f$ J_q\f$ forcing before the AD run.
  JqTermTLAD<MODEL> * initializeJqAD(const Increment_ &) const override;

/// Multiply by \f$ B\f$ and \f$ B^{-1}\f$.
  void Bmult(const Increment_ &, Increment_ &) const override;
  void Bminv(const Increment_ &, Increment_ &) const override;

/// Randomize
  void randomize(Increment_ &) const override;

/// Create new increment (set to 0).
  Increment_ * newStateIncrement() const override;

 private:
  std::unique_ptr<ModelSpaceCovarianceBase<MODEL> > B_;
  const State_ & xb_;
  const Variables ctlvars_;
  std::unique_ptr<const Geometry_> resol_;
  const eckit::LocalConfiguration conf_;
  const eckit::mpi::Comm & commTime_;
  const bool first_;
};

// =============================================================================

//  Generalized Jb Term of Cost Function
// -----------------------------------------------------------------------------

template<typename MODEL>
CostJbJq<MODEL>::CostJbJq(const eckit::Configuration & config, const eckit::mpi::Comm & comm,
                          const Geometry_ & resolouter, const Variables & ctlvars,
                          const State_ & xb)
  : B_(), xb_(xb), ctlvars_(ctlvars), resol_(), conf_(config), commTime_(comm),
    first_(comm.rank() == 0)
{
  Log::trace() << "CostJbJq contructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbJq<MODEL>::linearize(const State_ & fg, const Geometry_ & lowres) {
  Log::trace() << "CostJbJq::linearize start" << std::endl;
  resol_.reset(new Geometry_(lowres));
  const eckit::LocalConfiguration covConf(conf_, "background error");

  std::vector<eckit::LocalConfiguration> confs;
  covConf.get("covariances", confs);
  ASSERT(confs.size() == lowres.timeComm().size());
  eckit::LocalConfiguration myconf = confs[lowres.timeComm().rank()];

  B_.reset(CovarianceFactory<MODEL>::create(lowres, ctlvars_, myconf, xb_, fg));
  Log::trace() << "CostJbJq::linearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbJq<MODEL>::computeIncrement(const State_ & xb, const State_ & fg, const State_ & mx,
                                       Increment_ & dx) const {
  Log::trace() << "CostJbJq::computeIncrement start" << std::endl;
  static int tag = 13579;
  size_t mytime = commTime_.rank();
  State_ mxim1(fg);

// Send values of M(x_i) at end of my subwindow to next subwindow
  if (mytime + 1 < commTime_.size()) {
    oops::mpi::send(commTime_, mx, mytime+1, tag);
  }

// Receive values at beginning of my subwindow from previous subwindow
  if (mytime > 0) {
    oops::mpi::receive(commTime_, mxim1, mytime-1, tag);
  } else {
    mxim1 = xb;
  }

// Compute x_i - M(x_{i-1})
  dx.diff(fg, mxim1);

  ++tag;
  Log::info() << "CostJbJq: x_i - M(x_{i-1})" << dx << std::endl;
  Log::trace() << "CostJbJq::computeIncrement done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbJq<MODEL>::addGradient(const Increment_ & dxFG, Increment_ & grad,
                                  Increment_ & gradJb) const {
  Log::trace() << "CostJbJq::addGradient start" << std::endl;
  if (first_) {
//  Jb from pre-computed gradient
    grad += gradJb;

  } else {
//  Compute and add Jq gradient Qi^{-1} ( x_i - M(x_{i-1}) )
    Increment_ gg(grad, false);
    B_->inverseMultiply(dxFG, gg);
    grad += gg;
  }
  Log::trace() << "CostJbJq::addGradient done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
JqTermTLAD<MODEL> * CostJbJq<MODEL>::initializeJqTLAD() const {
  Log::trace() << "CostJbJq::initializeJqTLAD" << std::endl;
  return new JqTermTLAD<MODEL>(commTime_);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
JqTermTLAD<MODEL> * CostJbJq<MODEL>::initializeJqTL() const {
  Log::trace() << "CostJbJq::initializeJqTL start" << std::endl;
  JqTermTLAD<MODEL> * jqtl = new JqTermTLAD<MODEL>(commTime_);
  Log::trace() << "CostJbJq::initializeJqTL done" << std::endl;
  return jqtl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
JqTermTLAD<MODEL> * CostJbJq<MODEL>::initializeJqAD(const Increment_ & dx) const {
  Log::trace() << "CostJbJq::initializeJqAD start" << std::endl;
  JqTermTLAD<MODEL> * jqad = new JqTermTLAD<MODEL>(commTime_);
  jqad->setupAD(dx);
  Log::trace() << "CostJbJq::initializeJqAD done" << std::endl;
  return jqad;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbJq<MODEL>::Bmult(const Increment_ & dxin, Increment_ & dxout) const {
  Log::trace() << "CostJbJq::Bmult start" << std::endl;
  B_->multiply(dxin, dxout);
  Log::trace() << "CostJbJq::Bmult done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbJq<MODEL>::Bminv(const Increment_ & dxin, Increment_ & dxout) const {
  Log::trace() << "CostJbJq::Bminv start" << std::endl;
  B_->inverseMultiply(dxin, dxout);
  Log::trace() << "CostJbJq::Bminv done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbJq<MODEL>::randomize(Increment_ & dx) const {
  Log::trace() << "CostJbJq::randomize start" << std::endl;
  B_->randomize(dx);
  Log::trace() << "CostJbJq::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> * CostJbJq<MODEL>::newStateIncrement() const {
  Log::trace() << "CostJbJq::newStateIncrement start" << std::endl;
  Increment_ * incr = new Increment_(*resol_, ctlvars_, xb_.validTime());
  Log::trace() << "CostJbJq::newStateIncrement done" << std::endl;
  return incr;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJBJQ_H_
