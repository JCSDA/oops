/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2021-2023 UCAR
 * (C) Crown Copyright 2023, the Met Office.
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

#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostJbModelAux.h"
#include "oops/assimilation/CostJbObsAux.h"
#include "oops/assimilation/CostJbState.h"
#include "oops/base/Geometry.h"
#include "oops/base/ObsAuxCovariances.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/State.h"
#include "oops/interface/ModelAuxCovariance.h"
#include "oops/util/Logger.h"

namespace oops {
  template<typename MODEL> class JqTerm;
  template<typename MODEL> class JqTermTLAD;

// -----------------------------------------------------------------------------

/// Total Jb cost function for all components of the control variable.

template<typename MODEL, typename OBS> class CostJbTotal {
  typedef ControlIncrement<MODEL, OBS>  CtrlInc_;
  typedef ControlVariable<MODEL, OBS>   CtrlVar_;
  typedef State<MODEL>                  State_;
  typedef CostJbState<MODEL, OBS>       JbState_;
  typedef CostJbModelAux<MODEL, OBS>    JbModelAux_;
  typedef CostJbObsAux<MODEL, OBS>      JbObsAux_;
  typedef JqTerm<MODEL>                 JqTerm_;
  typedef JqTermTLAD<MODEL>             JqTermTLAD_;
  typedef Geometry<MODEL>               Geometry_;
  typedef ModelAuxCovariance<MODEL>     ModelAuxCovariance_;
  typedef ObsSpaces<OBS>                ObsSpaces_;
  typedef PostProcessor<State_>         PostProc_;
  typedef PostProcessorTLAD<MODEL>      PostProcTLAD_;
  typedef ObsAuxCovariances<OBS>        ObsAuxCovars_;

 public:
/// Construct \f$ J_b\f$.
  CostJbTotal(JbState_ *, const eckit::Configuration &, const Geometry_ &, const ObsSpaces_ &);

/// Destructor
  ~CostJbTotal();

/// Initialize before nonlinear model integration.
  void setPostProc(const CtrlVar_ &, const eckit::Configuration &, PostProc_ &);
  void setPostProcTraj(const CtrlVar_ &, const eckit::Configuration &,
                       const Geometry_ &, PostProcTLAD_ &);

/// Finalize computation after nonlinear model integration.
  double computeCost();
  void computeCostTraj();

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

/// Reset Jb gradient at first guess to zero
  void zeroGradientFG() {dxFG_->zero();}

/// Add Jb gradient at first guess.
  void addGradientFG(CtrlInc_ &) const;
  void addGradientFG(CtrlInc_ &, CtrlInc_ &) const;

/// Return background.
  const CtrlVar_ & getBackground() const {return xb_;}
  CtrlVar_ & getBackground() {return xb_;}

/// Return first guess \f$ x_0-x_b\f$.
  const CtrlInc_ & getFirstGuess() const {return *dxFG_;}

/// Jb terms for ControlIncrement constructor.
  const Geometry_ & resolution() const;
  const JbState_ & jbState() const {return *jb_;}
  const JbModelAux_ & jbModBias() const {return jbModBias_;}
  const JbObsAux_ & jbObsBias() const {return jbObsBias_;}
  const util::DateTime & windowBegin() const {return timeWindow_.start();}
  const util::DateTime & windowEnd()   const {return timeWindow_.end();}

 private:
  double evaluate(const CtrlInc_ &) const;

  std::unique_ptr<JbState_> jb_;
  JbModelAux_ jbModBias_;
  JbObsAux_  jbObsBias_;
  CtrlVar_ xb_;

/// Inner loop resolution
  const Geometry_ * resol_;
  const util::TimeWindow timeWindow_;
  eckit::LocalConfiguration innerConf_;

/// First guess increment \f$x_0-x_b\f$ or more generally \f$ x_i-M(x_{i-1})\f$.
  std::unique_ptr<CtrlInc_> dxFG_;
  const CtrlVar_ * xx_;
  const CtrlVar_ * traj_;

  bool jbEvaluation_;
  std::shared_ptr<JqTermTLAD_> jqtraj_;
  mutable std::shared_ptr<JqTermTLAD_> jqtl_;
  mutable std::shared_ptr<JqTermTLAD_> jqad_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostJbTotal<MODEL, OBS>::CostJbTotal(JbState_ * jb,
                                     const eckit::Configuration & conf,
                                     const Geometry_ & resol, const ObsSpaces_ & odb)
  : jb_(jb), jbModBias_(conf, resol),
    jbObsBias_(odb, conf.getSubConfiguration("observations.observers")),
    xb_(jb_->background(), jbModBias_.background(), jbObsBias_.background()),
    resol_(nullptr), timeWindow_(eckit::LocalConfiguration(conf, "time window")),
    innerConf_(), dxFG_(), xx_(), traj_(), jqtraj_()
{
  jbEvaluation_ = conf.getBool("jb evaluation", true);

  Log::trace() << "CostJbTotal contructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
CostJbTotal<MODEL, OBS>::~CostJbTotal() {
  Log::trace() << "CostJbTotal::~CostJbTotal start" << std::endl;
  jbObsBias_.covariance().write(jbObsBias_.covariance().config());
  Log::trace() << "CostJbTotal::~CostJbTotal done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
const Geometry<MODEL> & CostJbTotal<MODEL, OBS>::resolution() const {
  ASSERT(resol_);
  return *resol_;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbTotal<MODEL, OBS>::setPostProc(const CtrlVar_ & xx,
                                          const eckit::Configuration & config,
                                          PostProc_ & pp) {
  Log::trace() << "CostJbTotal::setPostProc start" << std::endl;
  xx_ = &xx;
  if (config.has("control pert")) jb_->setTime(xx);
  jb_->setPostProc(pp);
  Log::trace() << "CostJbTotal::setPostProc done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
double CostJbTotal<MODEL, OBS>::computeCost() {
  Log::trace() << "CostJbTotal::computeCost start" << std::endl;
  ASSERT(resol_);  // check that the B matrix has been setup (by setPostProcTraj), not great
  ASSERT(xx_);

  CtrlInc_ dx(*this);

// Compute x_0 - x_b for Jb (and Jq if present)
  jb_->computeIncrement(xb_, *xx_, jb_->getJq(), dx);

// Model and Obs biases
  dx.modVar().diff(xx_->modVar(), xb_.modVar());
  dx.obsVar().diff(xx_->obsVar(), xb_.obsVar());

// Print increment
  Log::info() << "CostJb: FG-BG" << dx << std::endl;

// Compute Jb value
  double zjb = 0.0;
  if (jbEvaluation_) zjb = this->evaluate(dx);

  xx_ = nullptr;
  Log::trace() << "CostJbTotal::computeCost done" << std::endl;
  return zjb;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbTotal<MODEL, OBS>::setPostProcTraj(const CtrlVar_ & fg,
                                              const eckit::Configuration & inner,
                                              const Geometry_ & resol, PostProcTLAD_ & pptraj) {
  Log::trace() << "CostJbTotal::setPostProcTraj start" << std::endl;
  traj_ = &fg;
  innerConf_ = eckit::LocalConfiguration(inner);
  resol_ = &resol;
// Trajectory for model error term
  jqtraj_.reset(jb_->initializeJqTLAD());
  pptraj.enrollProcessor(jqtraj_);
  jbModBias_.setPostProcTraj(*traj_, innerConf_, *resol_, pptraj);
  jbObsBias_.setPostProcTraj(*traj_, innerConf_, *resol_, pptraj);
  Log::trace() << "CostJbTotal::setPostProcTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbTotal<MODEL, OBS>::computeCostTraj() {
  Log::trace() << "CostJbTotal::computeCostTraj start" << std::endl;
  ASSERT(resol_);
  ASSERT(traj_);

// Linearize model-related terms and setup B (obs term done in computeCostTraj)
  jb_->linearize(xb_, *traj_, *resol_);
  jbModBias_.computeCostTraj();

// Linearize obs bias term
  jbObsBias_.computeCostTraj();

// Compute and save first guess increment.
  dxFG_.reset(new CtrlInc_(*this));

// Compute x_0 - x_b for Jb (and Jq if present)
  std::shared_ptr<JqTerm_> jq;
  if (jqtraj_) jq = jqtraj_->getJq();
  jb_->computeIncrement(xb_, *traj_, jq, *dxFG_);

// Model and Obs biases
  dxFG_->modVar().diff(traj_->modVar(), xb_.modVar());
  dxFG_->obsVar().diff(traj_->obsVar(), xb_.obsVar());

// Print increment
  Log::info() << "CostJb: FG-BG" << *dxFG_ << std::endl;
  jqtraj_.reset();
  traj_ = nullptr;
  Log::trace() << "CostJbTotal::computeCostTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
double CostJbTotal<MODEL, OBS>::evaluate(const CtrlInc_ & dx) const {
  Log::trace() << "CostJbTotal::evaluate start" << std::endl;
  ASSERT(resol_);
  CtrlInc_ gg(*this);
  this->multiplyBinv(dx, gg);

  double zjb = 0.0;
  double zz = 0.5 * dot_product(dx.states(), gg.states());
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
  jb_->addGradient(*dxFG_, grad, gradJb);
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
  jqad_.reset(jb_->initializeJqAD(dx));
  bgns += dx;
  ppad.enrollProcessor(jqad_);
  Log::trace() << "CostJbTotal::initializeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbTotal<MODEL, OBS>::finalizeAD() const {
  Log::trace() << "CostJbTotal::finalizeAD start" << std::endl;
  jqad_.reset();
  Log::trace() << "CostJbTotal::finalizeAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbTotal<MODEL, OBS>::multiplyB(const CtrlInc_ & dxin, CtrlInc_ & dxout) const {
  Log::trace() << "CostJbTotal::multiplyB start" << std::endl;
  jb_->Bmult(dxin, dxout);
  jbModBias_.multiplyCovar(dxin, dxout);
  jbObsBias_.multiplyCovar(dxin, dxout);
  Log::trace() << "CostJbTotal::multiplyB done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbTotal<MODEL, OBS>::multiplyBinv(const CtrlInc_ & dxin, CtrlInc_ & dxout) const {
  Log::trace() << "CostJbTotal::multiplyBinv start" << std::endl;
  jb_->Bminv(dxin, dxout);
  jbModBias_.multiplyCoInv(dxin, dxout);
  jbObsBias_.multiplyCoInv(dxin, dxout);
  Log::trace() << "CostJbTotal::multiplyBinv done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJbTotal<MODEL, OBS>::randomize(CtrlInc_ & dx) const {
  Log::trace() << "CostJbTotal::randomize start" << std::endl;
  jb_->randomize(dx);
  jbModBias_.randomize(dx);
  jbObsBias_.randomize(dx);
  Log::trace() << "CostJbTotal::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJBTOTAL_H_
