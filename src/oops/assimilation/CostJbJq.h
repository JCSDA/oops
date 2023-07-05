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

#ifndef OOPS_ASSIMILATION_COSTJBJQ_H_
#define OOPS_ASSIMILATION_COSTJBJQ_H_

#include <memory>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/CostJbState.h"
#include "oops/assimilation/JqTerm.h"
#include "oops/assimilation/JqTermTLAD.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment4D.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/State.h"
#include "oops/base/State4D.h"
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
  typedef Increment4D<MODEL>         Increment_;
  typedef State4D<MODEL>             State_;
  typedef JqTerm<MODEL>              JqTerm_;
  typedef JqTermTLAD<MODEL>          JqTLAD_;

 public:
/// Construct \f$ J_b\f$.
  CostJbJq(const std::vector<util::DateTime> &,
           const eckit::Configuration &, const eckit::mpi::Comm &,
           const Geometry_ &, const Variables &);

/// Destructor
  virtual ~CostJbJq() {}

  void setPostProc(PostProcessor<State<MODEL>> &) override;
  std::shared_ptr<JqTerm_> getJq() override {return jq_;}

/// Get increment from state (usually first guess).
  void computeIncrement(const State_ &, const State_ &, const std::shared_ptr<JqTerm_>,
                        Increment_ &) const override;

/// Linearize before the linear computations.
  void linearize(const State_ &, const State_ &, const Geometry_ &) override;

/// Add Jb gradient.
  void addGradient(const Increment_ &, Increment_ &, Increment_ &) const override;

/// Finalize \f$ J_q\f$ after the model run.
  JqTLAD_ * initializeJqTLAD() const override;

/// Finalize \f$ J_q\f$ after the TL run.
  JqTLAD_ * initializeJqTL() const override;

/// Initialize \f$ J_q\f$ forcing before the AD run.
  JqTLAD_ * initializeJqAD(const Increment_ &) const override;

/// Multiply by \f$ B\f$ and \f$ B^{-1}\f$.
  void Bmult(const Increment_ &, Increment_ &) const override;
  void Bminv(const Increment_ &, Increment_ &) const override;

/// Randomize
  void randomize(Increment_ &) const override;

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
  const eckit::LocalConfiguration conf_;
  const eckit::mpi::Comm & commTime_;
  std::shared_ptr<JqTerm_> jq_;
  std::vector<util::DateTime> times_;
};

// -----------------------------------------------------------------------------
//  Generalized Jb Term of Cost Function
// -----------------------------------------------------------------------------

template<typename MODEL>
CostJbJq<MODEL>::CostJbJq(const std::vector<util::DateTime> & times,
                          const eckit::Configuration & config, const eckit::mpi::Comm & comm,
                          const Geometry_ & geom, const Variables & ctlvars)
  : B_(), bg_(), ctlvars_(ctlvars), resol_(), conf_(config), commTime_(comm), jq_(), times_(times)
{
  Log::trace() << "CostJbJq::CostJbJq start" << std::endl;
  bg_.reset(new State_(times_, commTime_)),
  bg_->read(geom, eckit::LocalConfiguration(config, "background"));
  ASSERT(bg_->is_4d());
  Log::trace() << "CostJbJq::CostJbJq done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbJq<MODEL>::setPostProc(PostProcessor<State<MODEL>> & pp) {
  jq_.reset(new JqTerm_());
  pp.enrollProcessor(jq_);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbJq<MODEL>::linearize(const State_ & xb, const State_ & fg, const Geometry_ & lowres) {
  Log::trace() << "CostJbJq::linearize start" << std::endl;
  resol_ = &lowres;
  const eckit::LocalConfiguration covConf(conf_, "background error");

  std::vector<eckit::LocalConfiguration> confs;
  covConf.get("covariances", confs);
  ASSERT(confs.size() == times_.size());
  eckit::LocalConfiguration myconf = confs[commTime_.rank()];

  B_.reset(CovarianceFactory<MODEL>::create(lowres, ctlvars_, myconf, xb[0], fg[0]));
  Log::trace() << "CostJbJq::linearize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbJq<MODEL>::computeIncrement(const State_ & xb, const State_ & fg,
                                       const std::shared_ptr<JqTerm_> jq, Increment_ & dx) const {
  Log::trace() << "CostJbJq::computeIncrement start" << std::endl;
  ASSERT(jq);
  static int tag = 13579;
  size_t mytime = commTime_.rank();
  State_ mxim1(fg);

// Send values of M(x_i) at end of my subwindow to next subwindow
  if (mytime + 1 < commTime_.size()) {
    oops::mpi::send(commTime_, jq->getMofX(), mytime+1, tag);
  }

// Receive values at beginning of my subwindow from previous subwindow
  if (mytime > 0) {
    oops::mpi::receive(commTime_, mxim1[0], mytime-1, tag);
  } else {
    mxim1[0] = xb[0];
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
  if (commTime_.rank() == 0) {
//  Jb from pre-computed gradient
    grad += gradJb;

  } else {
//  Compute and add Jq gradient Qi^{-1} ( x_i - M(x_{i-1}) )
    Increment_ gg(grad, false);
    B_->inverseMultiply(dxFG[0], gg[0]);
    grad += gg;
  }
  Log::trace() << "CostJbJq::addGradient done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
JqTermTLAD<MODEL> * CostJbJq<MODEL>::initializeJqTLAD() const {
  Log::trace() << "CostJbJq::initializeJqTLAD" << std::endl;
  return new JqTLAD_(commTime_);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
JqTermTLAD<MODEL> * CostJbJq<MODEL>::initializeJqTL() const {
  Log::trace() << "CostJbJq::initializeJqTL start" << std::endl;
  JqTLAD_ * jqtl = new JqTLAD_(commTime_);
  Log::trace() << "CostJbJq::initializeJqTL done" << std::endl;
  return jqtl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
JqTermTLAD<MODEL> * CostJbJq<MODEL>::initializeJqAD(const Increment_ & dx) const {
  Log::trace() << "CostJbJq::initializeJqAD start" << std::endl;
  JqTLAD_ * jqad = new JqTLAD_(commTime_);
  jqad->setupAD(dx[0]);
  Log::trace() << "CostJbJq::initializeJqAD done" << std::endl;
  return jqad;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbJq<MODEL>::Bmult(const Increment_ & dxin, Increment_ & dxout) const {
  Log::trace() << "CostJbJq::Bmult start" << std::endl;
  B_->multiply(dxin[0], dxout[0]);
  Log::trace() << "CostJbJq::Bmult done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbJq<MODEL>::Bminv(const Increment_ & dxin, Increment_ & dxout) const {
  Log::trace() << "CostJbJq::Bminv start" << std::endl;
  B_->inverseMultiply(dxin[0], dxout[0]);
  Log::trace() << "CostJbJq::Bminv done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbJq<MODEL>::randomize(Increment_ & dx) const {
  Log::trace() << "CostJbJq::randomize start" << std::endl;
  B_->randomize(dx[0]);
  Log::trace() << "CostJbJq::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJBJQ_H_
