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
#include "oops/base/Geometry.h"
#include "oops/base/ObsAuxCovariances.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/State.h"
#include "oops/interface/ModelAuxCovariance.h"
#include "oops/util/Logger.h"

namespace oops {
  template<typename MODEL> class JqTermTLAD;

// -----------------------------------------------------------------------------

/// Total Jb cost function for all components of the control variable.

template<typename MODEL, typename OBS> class CostJbTotal {
  typedef ControlIncrement<MODEL, OBS>  CtrlInc_;
  typedef ControlVariable<MODEL, OBS>   CtrlVar_;
  typedef State<MODEL>                  State_;
  typedef CostJbState<MODEL>            JbState_;
  typedef JqTermTLAD<MODEL>             JqTermTLAD_;
  typedef Geometry<MODEL>               Geometry_;
  typedef ModelAuxCovariance<MODEL>     ModelAuxCovar_;
  typedef ObsAuxCovariances<OBS>        ObsAuxCovars_;
  typedef ObsSpaces<OBS>                ObsSpaces_;
  typedef PostProcessorTLAD<MODEL>      PostProcTLAD_;

 public:
/// Construct \f$ J_b\f$.
  CostJbTotal(const CtrlVar_ &, JbState_ *, const eckit::Configuration &,
              const Geometry_ &, const ObsSpaces_ & odb);

/// Destructor
  ~CostJbTotal() {}

/// Initialize before nonlinear model integration.
  void initialize(const CtrlVar_ &) const;
  void initializeTraj(const CtrlVar_ &, const Geometry_ &,
                      const eckit::Configuration &, PostProcTLAD_ &);

/// Finalize computation after nonlinear model integration.
  double finalize(const CtrlVar_ &) const;
  void finalizeTraj();

/// Initialize before starting the TL run.
  void initializeTL(PostProcTLAD_ &) const;
  void finalizeTL(const CtrlInc_ &, CtrlInc_ &) const;

/// Initialize before starting the AD run.
  void initializeAD(CtrlInc_ &, const CtrlInc_ &, PostProcTLAD_ &) const;
  void finalizeAD() const;

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
  const util::DateTime & windowBegin() const {return windowBegin_;}
  const util::DateTime & windowEnd()   const {return windowEnd_;}

 private:
  double evaluate(const CtrlInc_ &) const;

  const CtrlVar_ & xb_;
  std::unique_ptr<JbState_> jb_;
  ModelAuxCovar_ jbModBias_;
  ObsAuxCovars_  jbObsBias_;

  const CtrlVar_ mutable * fg_;

/// First guess increment \f$x_0-x_b\f$ or more generally \f$ x_i-M(x_{i-1})\f$.
  std::unique_ptr<CtrlInc_> dxFG_;

/// Inner loop resolution
  std::unique_ptr<Geometry_> resol_;
  const util::DateTime windowBegin_;
  const util::DateTime windowEnd_;

  bool jbEvaluation_;
  std::shared_ptr<JqTermTLAD_> jqtraj_;
  mutable std::shared_ptr<JqTermTLAD_> jqtl_;
  mutable std::shared_ptr<JqTermTLAD_> jqad_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostJbTotal<MODEL, OBS>::CostJbTotal(const CtrlVar_ & xb, JbState_ * jb,
                                     const eckit::Configuration & conf,
                                     const Geometry_ & resol, const ObsSpaces_ & odb)
  : xb_(xb), jb_(jb),
    jbModBias_(conf.getSubConfiguration("model aux error"), resol),
    jbObsBias_(odb, conf.getSubConfiguration("observations")), dxFG_(), resol_(),
    windowBegin_(conf.getString("window begin")),
    windowEnd_(windowBegin_ + util::Duration(conf.getString("window length"))),
    jqtraj_()
{
  jbEvaluation_ = conf.getBool("jb evaluation", true);
  Log::trace() << "CostJbTotal contructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbTotal<MODEL, OBS>::initialize(const CtrlVar_ & fg) const {
  Log::trace() << "CostJbTotal::initialize start" << std::endl;
  fg_ = &fg;
  Log::trace() << "CostJbTotal::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
double CostJbTotal<MODEL, OBS>::finalize(const CtrlVar_ & mx) const {
  Log::trace() << "CostJbTotal::finalize start" << std::endl;
  ASSERT(fg_);
  CtrlInc_ dx(*this);

// Compute x_0 - x_b for Jb (and Jq is present)
  jb_->computeIncrement(xb_.state(), fg_->state(), mx.state(), dx.state());

// Model and Obs biases
  dx.modVar().diff(fg_->modVar(), xb_.modVar());
  dx.obsVar().diff(fg_->obsVar(), xb_.obsVar());

// Print increment
  Log::info() << "CostJb: FG-BG" << dx << std::endl;

// Compute Jb value
  double zjb = 0.0;
  if (jbEvaluation_) zjb = this->evaluate(dx);
  Log::trace() << "CostJbTotal::finalize done" << std::endl;
  return zjb;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbTotal<MODEL, OBS>::initializeTraj(const CtrlVar_ & fg, const Geometry_ & resol,
                                             const eckit::Configuration & inner,
                                             PostProcTLAD_ & pptraj) {
  Log::trace() << "CostJbTotal::initializeTraj start" << std::endl;
  fg_ = &fg;
  resol_.reset(new Geometry_(resol));
// Linearize terms
  jb_->linearize(fg.state(), *resol_);
  jbModBias_.linearize(fg.modVar(), *resol_);
  jbObsBias_.linearize(fg.obsVar(), inner);

  jqtraj_.reset(jb_->initializeJqTLAD());
  pptraj.enrollProcessor(jqtraj_);

  Log::trace() << "CostJbTotal::initializeTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbTotal<MODEL, OBS>::finalizeTraj() {
  Log::trace() << "CostJbTotal::finalizeTraj start" << std::endl;
  ASSERT(fg_);
// Compute and save first guess increment.
  dxFG_.reset(new CtrlInc_(*this));

// Compute x_0 - x_b for Jb (and Jq is present)
  const State_ * mx = &fg_->state();
  if (jqtraj_) mx = &jqtraj_->getMxi();
  jb_->computeIncrement(xb_.state(), fg_->state(), *mx, dxFG_->state());

// Model and Obs biases
  dxFG_->modVar().diff(fg_->modVar(), xb_.modVar());
  dxFG_->obsVar().diff(fg_->obsVar(), xb_.obsVar());

// Print increment
  Log::info() << "CostJb: FG-BG" << *dxFG_ << std::endl;
  jqtraj_.reset();
  Log::trace() << "CostJbTotal::finalizeTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
double CostJbTotal<MODEL, OBS>::evaluate(const CtrlInc_ & dx) const {
  Log::trace() << "CostJbTotal::evaluate start" << std::endl;
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

  Log::trace() << "CostJbTotal::evaluate done" << std::endl;
  return zjb;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbTotal<MODEL, OBS>::addGradientFG(CtrlInc_ & grad) const {
  Log::trace() << "CostJbTotal::addGradientFG 1 start" << std::endl;
  CtrlInc_ gg(grad, false);
  this->multiplyBinv(*dxFG_, gg);
  grad += gg;
  Log::trace() << "CostJbTotal::addGradientFG 1 done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbTotal<MODEL, OBS>::addGradientFG(CtrlInc_ & grad, CtrlInc_ & gradJb) const {
  Log::trace() << "CostJbTotal::addGradientFG 2 start" << std::endl;
  jb_->addGradient(dxFG_->state(), grad.state(), gradJb.state());
  grad.modVar() += gradJb.modVar();
  grad.obsVar() += gradJb.obsVar();
  Log::trace() << "CostJbTotal::addGradientFG 2 done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbTotal<MODEL, OBS>::initializeTL(PostProcTLAD_ & pptl) const {
  Log::trace() << "CostJbTotal::initializeTL start" << std::endl;
  jqtl_.reset(jb_->initializeJqTL());
  pptl.enrollProcessor(jqtl_);
  Log::trace() << "CostJbTotal::initializeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbTotal<MODEL, OBS>::finalizeTL(const CtrlInc_ & bgns, CtrlInc_ & dx) const {
  Log::trace() << "CostJbTotal::finalizeTL start" << std::endl;
  dx = bgns;
  if (jqtl_) jqtl_->computeModelErrorTL(dx.state());
  jqtl_.reset();
  Log::trace() << "CostJbTotal::finalizeTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbTotal<MODEL, OBS>::initializeAD(CtrlInc_ & bgns, const CtrlInc_ & dx,
                                           PostProcTLAD_ & ppad) const {
  Log::trace() << "CostJbTotal::initializeAD start" << std::endl;
  jqad_.reset(jb_->initializeJqAD(dx.state()));
  bgns += dx;
  ppad.enrollProcessor(jqad_);
  Log::trace() << "CostJbTotal::initializeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbTotal<MODEL, OBS>::finalizeAD() const {
  Log::trace() << "CostJbTotal::finalizeAD start" << std::endl;
  if (jqad_) jqad_->clear();
  jqad_.reset();
  Log::trace() << "CostJbTotal::finalizeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbTotal<MODEL, OBS>::multiplyB(const CtrlInc_ & dxin, CtrlInc_ & dxout) const {
  Log::trace() << "CostJbTotal::multiplyB start" << std::endl;
  jb_->Bmult(dxin.state(), dxout.state());
  jbModBias_.multiply(dxin.modVar(), dxout.modVar());
  jbObsBias_.multiply(dxin.obsVar(), dxout.obsVar());
  Log::trace() << "CostJbTotal::multiplyB done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbTotal<MODEL, OBS>::multiplyBinv(const CtrlInc_ & dxin, CtrlInc_ & dxout) const {
  Log::trace() << "CostJbTotal::multiplyBinv start" << std::endl;
  jb_->Bminv(dxin.state(), dxout.state());
  jbModBias_.inverseMultiply(dxin.modVar(), dxout.modVar());
  jbObsBias_.inverseMultiply(dxin.obsVar(), dxout.obsVar());
  Log::trace() << "CostJbTotal::multiplyBinv done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbTotal<MODEL, OBS>::randomize(CtrlInc_ & dx) const {
  Log::trace() << "CostJbTotal::randomize start" << std::endl;
  jb_->randomize(dx.state());
  jbModBias_.randomize(dx.modVar());
  jbObsBias_.randomize(dx.obsVar());
  Log::trace() << "CostJbTotal::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJBTOTAL_H_
