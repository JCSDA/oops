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

#include<limits>

#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "util/Logger.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostJbState.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/ModelAuxCovariance.h"
#include "oops/interface/ObsAuxCovariance.h"

namespace oops {
  template<typename MODEL> class JqTerm;
  template<typename MODEL> class JqTermTL;
  template<typename MODEL> class JqTermAD;

// -----------------------------------------------------------------------------

/// Total Jb cost function for all components of the control variable.

template<typename MODEL> class CostJbTotal {
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef ControlVariable<MODEL>     CtrlVar_;
  typedef CostJbState<MODEL>         JbState_;
  typedef JqTerm<MODEL>              JqTerm_;
  typedef JqTermAD<MODEL>            JqTermAD_;
  typedef JqTermTL<MODEL>            JqTermTL_;
  typedef Geometry<MODEL>            Geometry_;
  typedef ModelAuxCovariance<MODEL>  ModelAuxCovar_;
  typedef ObsAuxCovariance<MODEL>    ObsAuxCovar_;

 public:
/// Construct \f$ J_b\f$.
  CostJbTotal(const CtrlVar_ &, JbState_ *, const eckit::Configuration &, const Geometry_ &);

/// Destructor
  ~CostJbTotal() {}

/// Initialize before nonlinear model integration.
  JqTerm_ * initialize(const CtrlVar_ &) const;
  JqTerm_ * initializeTraj(const CtrlVar_ &, const Geometry_ &);

/// Finalize computation after nonlinear model integration.
  double finalize(JqTerm_ *) const;
  double finalizeTraj(JqTerm_ *);

/// Initialize before starting the TL run.
  JqTermTL_ * initializeTL() const;
  void finalizeTL(JqTermTL_ *, const CtrlInc_ &, CtrlInc_ &) const;

/// Initialize before starting the AD run.
  JqTermAD_ * initializeAD(CtrlInc_ &, const CtrlInc_ &) const;
  void finalizeAD(JqTermAD_ *) const;

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
  const ObsAuxCovar_ & jbObsBias() const {return jbObsBias_;}

 private:
  double evaluate(const CtrlInc_ &) const;

  const CtrlVar_ & xb_;
  boost::scoped_ptr<JbState_> jb_;
  ModelAuxCovar_ jbModBias_;
  ObsAuxCovar_   jbObsBias_;

  const CtrlVar_ mutable * fg_;

/// First guess increment \f$x_0-x_b\f$ or more generally \f$ x_i-M(x_{i-1})\f$.
  boost::scoped_ptr<CtrlInc_> dxFG_;

/// Inner loop resolution
  boost::scoped_ptr<Geometry_> resol_;
};

// =============================================================================

template<typename MODEL>
CostJbTotal<MODEL>::CostJbTotal(const CtrlVar_ & xb, JbState_ * jb,
                                const eckit::Configuration & conf, const Geometry_ & resol)
  : xb_(xb), jb_(jb),
    jbModBias_(conf.getSubConfiguration("ModelBiasCovariance"), resol),
    jbObsBias_(conf.getSubConfiguration("ObsBiasCovariance")),
    dxFG_(0), resol_(0)
{
  Log::trace() << "CostJbTotal contructed." << std::endl;
  Log::debug() << "CostJbTotal background is:" << xb_ << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
JqTerm<MODEL> * CostJbTotal<MODEL>::initialize(const CtrlVar_ & fg) const {
  fg_ = &fg;
  Log::debug() << "CostJbTotal first guess is:" << *fg_ << std::endl;
  JqTerm_ * jqnl = jb_->initializeJq();
  return jqnl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double CostJbTotal<MODEL>::finalize(JqTerm_ * jqnl) const {
  CtrlInc_ dx(*this);
  Log::debug() << "CostJbTotal:finalize background is:" << xb_ << std::endl;
  Log::debug() << "CostJbTotal:finalize first guess is:" << *fg_ << std::endl;

// Compute x_0 - x_b for Jb
  jb_->computeIncrement(xb_.state(), fg_->state(), dx.state());

// Model and Obs biases
  dx.modVar().diff(fg_->modVar(), xb_.modVar());
  dx.obsVar().diff(fg_->obsVar(), xb_.obsVar());

  if (jqnl) jqnl->computeModelError(fg_->state(), dx.state());
  Log::info() << "CostJb: FG-BG" << dx << std::endl;

  fg_ = NULL;

// Compute Jb value
  double zjb = this->evaluate(dx);
  return zjb;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
JqTerm<MODEL> * CostJbTotal<MODEL>::initializeTraj(const CtrlVar_ & fg,
                                                   const Geometry_ & resol) {
  resol_.reset(new Geometry_(resol));
// Linearize terms
  jb_->linearize(fg.state(), *resol_);
  jbModBias_.linearize(fg.modVar(), *resol_);
  jbObsBias_.linearize(fg.obsVar());

  JqTerm_ * jqnl = this->initialize(fg);
  if (jqnl) jqnl->linearize();
  return jqnl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double CostJbTotal<MODEL>::finalizeTraj(JqTerm_ * jqnl) {
// Compute and save first guess increment.
  dxFG_.reset(new CtrlInc_(*this));

  Log::debug() << "CostJbTotal:finalizeTraj background is:" << xb_ << std::endl;
  Log::debug() << "CostJbTotal:finalizeTraj first guess is:" << *fg_ << std::endl;

// Compute x_0 - x_b for Jb
  jb_->computeIncrement(xb_.state(), fg_->state(), dxFG_->state());

// Model and Obs biases
  dxFG_->modVar().diff(fg_->modVar(), xb_.modVar());
  dxFG_->obsVar().diff(fg_->obsVar(), xb_.obsVar());

  if (jqnl) jqnl->computeModelError(fg_->state(), dxFG_->state());
  Log::info() << "CostJb: FG-BG" << *dxFG_ << std::endl;

  fg_ = NULL;

// Return Jb value
  double zjb = this->evaluate(*dxFG_);
  return zjb;
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
JqTermTL<MODEL> * CostJbTotal<MODEL>::initializeTL() const {
  JqTermTL_ * jqtl = jb_->initializeJqTL();
  return jqtl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbTotal<MODEL>::finalizeTL(JqTermTL_ * jqtl, const CtrlInc_ & bgns,
                                    CtrlInc_ & dx) const {
  dx = bgns;
  if (jqtl) jqtl->computeModelErrorTL(dx.state());
}

// -----------------------------------------------------------------------------

template<typename MODEL>
JqTermAD<MODEL> * CostJbTotal<MODEL>::initializeAD(CtrlInc_ & bgns,
                                                   const CtrlInc_ & dx) const {
  JqTermAD_ * jqad = jb_->initializeJqAD(dx.state());
  bgns += dx;
  return jqad;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJbTotal<MODEL>::finalizeAD(JqTermAD_ * jqad) const {
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

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJBTOTAL_H_
