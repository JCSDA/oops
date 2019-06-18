/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_COSTJBTOTAL_H_
#define OOPS_ASSIMILATION_COSTJBTOTAL_H_

#include <limits>
#include <memory>


#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostJbState.h"
#include "oops/base/ObsAuxCovariances.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/ModelAuxCovariance.h"
#include "oops/util/Logger.h"

namespace oops {
  template<typename MODEL> class JqTerm;
  template<typename MODEL> class JqTermTLAD;

// -----------------------------------------------------------------------------

/// Total Jb cost function for all components of the control variable.

template<typename MODEL> class CostJbTotal {
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef ControlVariable<MODEL>     CtrlVar_;
  typedef CostJbState<MODEL>         JbState_;
  typedef JqTerm<MODEL>              JqTerm_;
  typedef JqTermTLAD<MODEL>          JqTermTLAD_;
  typedef Geometry<MODEL>            Geometry_;
  typedef ModelAuxCovariance<MODEL>  ModelAuxCovar_;
  typedef ObsAuxCovariances<MODEL>   ObsAuxCovars_;

 public:
/// Construct \f$ J_b\f$.
  CostJbTotal(const CtrlVar_ &, JbState_ *, const eckit::Configuration &, const Geometry_ &);

/// Destructor
  ~CostJbTotal() {}

/// Initialize before nonlinear model integration.
  JqTerm_ * initialize(const CtrlVar_ &) const;
  JqTermTLAD_ * initializeTraj(const CtrlVar_ &, const Geometry_ &);

/// Finalize computation after nonlinear model integration.
  double finalize(JqTerm_ *) const;
  void finalizeTraj(JqTermTLAD_ *);

/// Initialize before starting the TL run.
  JqTermTLAD_ * initializeTL() const;
  void finalizeTL(JqTermTLAD_ *, const CtrlInc_ &, CtrlInc_ &) const;

/// Initialize before starting the AD run.
  JqTermTLAD_ * initializeAD(CtrlInc_ &, const CtrlInc_ &) const;
  void finalizeAD(JqTermTLAD_ *) const;

/// Multiply by covariance matrix and its inverse.
  void multiplyB(const CtrlInc_ &, CtrlInc_ &) const;
  void multiplyBinv(const CtrlInc_ &, CtrlInc_ &) const;

/// Randomize
  void randomize(CtrlInc_ &) const;

/// Add Jb gradient at first guess.
  void addGradientFG(CtrlInc_ &) const;
  void addGradientFG(CtrlInc_ &, CtrlInc_ &) const;

/// Return background.
  const CtrlVar_ & getBackground() const {return xb_;}

/// Return first guess \f$ x_0-x_b\f$.
  const CtrlInc_ & getFirstGuess() const {return *dxFG_;}

/// Jb terms for ControlIncrement constructor.
  const Geometry_ & resolution() const {return *resol_;}
  const JbState_ & jbState() const {return *jb_;}
  const ModelAuxCovar_ & jbModBias() const {return jbModBias_;}
  const ObsAuxCovars_ & jbObsBias() const {return jbObsBias_;}

 private:
  double evaluate(const CtrlInc_ &) const;

  const CtrlVar_ & xb_;
  std::unique_ptr<JbState_> jb_;
  ModelAuxCovar_ jbModBias_;
  ObsAuxCovars_   jbObsBias_;

  const CtrlVar_ mutable * fg_;

/// First guess increment \f$x_0-x_b\f$ or more generally \f$ x_i-M(x_{i-1})\f$.
  std::unique_ptr<CtrlInc_> dxFG_;

/// Inner loop resolution
  std::unique_ptr<Geometry_> resol_;
};

// =============================================================================

template<typename MODEL>
CostJbTotal<MODEL>::CostJbTotal(const CtrlVar_ & xb, JbState_ * jb,
                                const eckit::Configuration & conf, const Geometry_ & resol)
  : xb_(xb), jb_(jb),
    jbModBias_(conf.getSubConfiguration("Jb.ModelBiasCovariance"), resol),
    jbObsBias_(eckit::LocalConfiguration(conf, "Jo")),
    dxFG_(), resol_()
{
  Log::trace() << "CostJbTotal contructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
JqTerm<MODEL> * CostJbTotal<MODEL>::initialize(const CtrlVar_ & fg) const {
  fg_ = &fg;
  JqTerm_ * jqnl = jb_->initializeJq();
  return jqnl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double CostJbTotal<MODEL>::finalize(JqTerm_ * jqnl) const {
  ASSERT(fg_);
  CtrlInc_ dx(*this);

// Compute x_0 - x_b for Jb
  jb_->computeIncrement(xb_.state(), fg_->state(), dx.state());

// Model and Obs biases
  dx.modVar().diff(fg_->modVar(), xb_.modVar());
  dx.obsVar().diff(fg_->obsVar(), xb_.obsVar());

  if (jqnl) jqnl->computeModelError(fg_->state(), dx.state());
  Log::info() << "CostJb: FG-BG" << dx << std::endl;

// Compute Jb value
  double zjb = this->evaluate(dx);
  return zjb;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
JqTermTLAD<MODEL> * CostJbTotal<MODEL>::initializeTraj(const CtrlVar_ & fg,
                                                       const Geometry_ & resol) {
  fg_ = &fg;
  resol_.reset(new Geometry_(resol));
// Linearize terms
  jb_->linearize(fg.state(), *resol_);
  jbModBias_.linearize(fg.modVar(), *resol_);
  jbObsBias_.linearize(fg.obsVar());

  JqTermTLAD_ * jqlin = jb_->initializeJqTLAD();

  return jqlin;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbTotal<MODEL>::finalizeTraj(JqTermTLAD_ * jqlin) {
  ASSERT(fg_);
// Compute and save first guess increment.
  dxFG_.reset(new CtrlInc_(*this));

// Compute x_0 - x_b for Jb
  jb_->computeIncrement(xb_.state(), fg_->state(), dxFG_->state());

// Model and Obs biases
  dxFG_->modVar().diff(fg_->modVar(), xb_.modVar());
  dxFG_->obsVar().diff(fg_->obsVar(), xb_.obsVar());

  if (jqlin) jqlin->computeModelErrorTraj(fg_->state(), dxFG_->state());
  Log::info() << "CostJb: FG-BG" << *dxFG_ << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double CostJbTotal<MODEL>::evaluate(const CtrlInc_ & dx) const {
  CtrlInc_ gg(*this);
  this->multiplyBinv(dx, gg);

  double zjb = 0.0;
  double zz = 0.5 * dot_product(dx.state(), gg.state());
  Log::info() << "CostJb   : Nonlinear Jb State = " << zz << std::endl;
  zjb += zz;
  zz = 0.5 * dot_product(dx.modVar(), gg.modVar());
  Log::info() << "CostJb   : Nonlinear Jb Model Aux = " << zz << std::endl;
  zjb += zz;
  zz = 0.5 * dot_product(dx.obsVar(), gg.obsVar());
  Log::info() << "CostJb   : Nonlinear Jb Obs Aux = " << zz << std::endl;
  zjb += zz;

  Log::info() << "CostJb   : Nonlinear Jb = " << zjb << std::endl;

// Get rid of very small values for test
  double ztest = zjb;
  if (zjb >= 0.0 && zjb <= std::numeric_limits<double>::epsilon()) ztest = 0.0;
  Log::test() << "CostJb   : Nonlinear Jb = " << ztest << std::endl;

  return zjb;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbTotal<MODEL>::addGradientFG(CtrlInc_ & grad) const {
  CtrlInc_ gg(grad, false);
  this->multiplyBinv(*dxFG_, gg);
  grad += gg;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbTotal<MODEL>::addGradientFG(CtrlInc_ & grad, CtrlInc_ & gradJb) const {
  jb_->addGradient(dxFG_->state(), grad.state(), gradJb.state());
  grad.modVar() += gradJb.modVar();
  grad.obsVar() += gradJb.obsVar();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
JqTermTLAD<MODEL> * CostJbTotal<MODEL>::initializeTL() const {
  JqTermTLAD_ * jqtl = jb_->initializeJqTL();
  return jqtl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbTotal<MODEL>::finalizeTL(JqTermTLAD_ * jqtl, const CtrlInc_ & bgns,
                                    CtrlInc_ & dx) const {
  dx = bgns;
  if (jqtl) jqtl->computeModelErrorTL(dx.state());
}

// -----------------------------------------------------------------------------

template<typename MODEL>
JqTermTLAD<MODEL> * CostJbTotal<MODEL>::initializeAD(CtrlInc_ & bgns,
                                                   const CtrlInc_ & dx) const {
  JqTermTLAD_ * jqad = jb_->initializeJqAD(dx.state());
  bgns += dx;
  return jqad;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbTotal<MODEL>::finalizeAD(JqTermTLAD_ * jqad) const {
  if (jqad) jqad->clear();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbTotal<MODEL>::multiplyB(const CtrlInc_ & dxin, CtrlInc_ & dxout) const {
  jb_->Bmult(dxin.state(), dxout.state());
  jbModBias_.multiply(dxin.modVar(), dxout.modVar());
  jbObsBias_.multiply(dxin.obsVar(), dxout.obsVar());
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbTotal<MODEL>::multiplyBinv(const CtrlInc_ & dxin, CtrlInc_ & dxout) const {
  jb_->Bminv(dxin.state(), dxout.state());
  jbModBias_.inverseMultiply(dxin.modVar(), dxout.modVar());
  jbObsBias_.inverseMultiply(dxin.obsVar(), dxout.obsVar());
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbTotal<MODEL>::randomize(CtrlInc_ & dx) const {
  jb_->randomize(dx.state());
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJBTOTAL_H_
