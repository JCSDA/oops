/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_COSTJCDFI_H_
#define OOPS_ASSIMILATION_COSTJCDFI_H_

#include <memory>
#include <utility>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/DolphChebyshev.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/base/WeightedDiff.h"
#include "oops/base/WeightedDiffTLAD.h"
#include "oops/base/WeightingFct.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Jc DFI Cost Function
/*!
 * Digital filter based constraint term for the cost function.
 */

template<typename MODEL, typename OBS> class CostJcDFI : public CostTermBase<MODEL, OBS> {
  typedef ControlIncrement<MODEL, OBS>  CtrlInc_;
  typedef ControlVariable<MODEL, OBS>   CtrlVar_;
  typedef Geometry<MODEL>               Geometry_;
  typedef Increment<MODEL>              Increment_;
  typedef State<MODEL>                  State_;
  typedef PostProcessor<State_>         PostProc_;
  typedef PostProcessorTLAD<MODEL>      PostProcTLAD_;

 public:
/// Construct \f$ J_c\f$.
  CostJcDFI(const eckit::Configuration &, const Geometry_ &, const util::DateTime &,
            const util::Duration &, const util::Duration & tstep = util::Duration(0));

/// Destructor
  virtual ~CostJcDFI() {}

/// Nonlinear Jc DFI computation
  void setPostProc(const CtrlVar_ &, const eckit::Configuration &, PostProc_ &) override;
  double computeCost() override;
  void printCostTestHack() override;

/// Linearization trajectory for Jc DFI computation
  void setPostProcTraj(const CtrlVar_ &, const eckit::Configuration &,
                       const Geometry_ &, PostProcTLAD_ &) override;
  void computeCostTraj() override;

/// TL Jc DFI computation
  void setPostProcTL(const CtrlInc_ &, PostProcTLAD_ &) const override;
  void computeCostTL(const CtrlInc_ &, GeneralizedDepartures &) const override;

/// Adjoint Jc DFI computation
  void computeCostAD(std::shared_ptr<const GeneralizedDepartures>,
                     CtrlInc_ &, PostProcTLAD_ &) const override;
  void setPostProcAD() const override {}

/// Multiply by \f$ C\f$ and \f$ C^{-1}\f$.
  std::unique_ptr<GeneralizedDepartures>
    multiplyCovar(const GeneralizedDepartures &) const override;
  std::unique_ptr<GeneralizedDepartures>
    multiplyCoInv(const GeneralizedDepartures &) const override;

/// Provide new increment.
  std::unique_ptr<GeneralizedDepartures> newDualVector() const override;

/// Gradient of \f$ J_c\f$ at first guess.
  std::unique_ptr<GeneralizedDepartures> newGradientFG() const override;

/// Reset trajectory.
  void resetLinearization() override;

 private:
  util::DateTime vt_;
  util::Duration span_;
  double alpha_;
  std::unique_ptr<WeightingFct> wfct_;
  std::unique_ptr<Increment_> gradFG_;
  const Geometry_ & resol_;
  const util::Duration tstep_;
  const Geometry_ * tlres_;
  util::Duration tlstep_;
  mutable std::shared_ptr<WeightedDiff<MODEL, Increment_, State_> > filter_;
  mutable std::shared_ptr<WeightedDiffTLAD<MODEL> > ftlad_;
  Variables vars_;
  double zhack_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostJcDFI<MODEL, OBS>::CostJcDFI(const eckit::Configuration & conf, const Geometry_ & resol,
                                 const util::DateTime & vt, const util::Duration & span,
                                 const util::Duration & tstep)
  : vt_(vt), span_(span), alpha_(0), wfct_(), gradFG_(),
    resol_(resol), tstep_(tstep), tlres_(), tlstep_(), filter_(),
    vars_(conf, "filtered variables"), zhack_(util::missingValue<double>())
{
  alpha_ = conf.getDouble("alpha");
  if (conf.has("ftime")) vt_ = util::DateTime(conf.getString("ftime"));
  if (conf.has("span")) span_ = util::Duration(conf.getString("span"));
//  wfct_.reset(WeightFactory::create(config)); YT
  wfct_.reset(new DolphChebyshev(conf));
  Log::trace() << "CostJcDFI created" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJcDFI<MODEL, OBS>::setPostProc(const CtrlVar_ &, const eckit::Configuration &,
                                        PostProc_ & pp) {
  filter_.reset(new WeightedDiff<MODEL, Increment_, State_>(vars_, vt_, span_,
                                                            tstep_, resol_, *wfct_));
  pp.enrollProcessor(filter_);
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
double CostJcDFI<MODEL, OBS>::computeCost() {
  double zz = 0.5 * alpha_;
  std::unique_ptr<Increment_> dx(filter_->releaseDiff());
  zz *= dot_product(*dx, *dx);
  Log::info() << "CostJcDFI: Nonlinear Jc = " << zz << std::endl;
// TEMPORARY HACK START
//  Log::test() << "CostJcDFI: Nonlinear Jc = " << zz << std::endl;
  zhack_ = zz;
// TEMPORARY HACK END
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJcDFI<MODEL, OBS>::printCostTestHack() {
  Log::test() << "CostJcDFI: Nonlinear Jc = " << zhack_ << std::endl;
  zhack_ = util::missingValue<double>();
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJcDFI<MODEL, OBS>::setPostProcTraj(const CtrlVar_ &, const eckit::Configuration & conf,
                                            const Geometry_ & tlres, PostProcTLAD_ & pptraj) {
  tlres_ = &tlres;
  tlstep_ = util::Duration(conf.getString("linear model.tstep", tstep_.toString()));
  ftlad_.reset(new WeightedDiffTLAD<MODEL>(vars_, vt_, span_, tstep_, *tlres_, *wfct_));
  pptraj.enrollProcessor(ftlad_);
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJcDFI<MODEL, OBS>::computeCostTraj() {
  gradFG_.reset(ftlad_->releaseDiff());
  *gradFG_ *= alpha_;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJcDFI<MODEL, OBS>::setPostProcTL(const CtrlInc_ &, PostProcTLAD_ & pptl) const {
  ftlad_->setupTL(*tlres_);
  pptl.enrollProcessor(ftlad_);
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJcDFI<MODEL, OBS>::computeCostTL(const CtrlInc_ & dx, GeneralizedDepartures & gdep) const {
  Log::trace() << "CostJcDFI::computeCostTL start" << std::endl;
  Increment_ & ydep = dynamic_cast<Increment_ &>(gdep);
  ftlad_->finalTL(ydep);
  Log::trace() << "CostJcDFI::computeCostTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJcDFI<MODEL, OBS>::computeCostAD(std::shared_ptr<const GeneralizedDepartures> pv,
                                          CtrlInc_ &, PostProcTLAD_ & ppad) const {
  Log::trace() << "CostJcDFI::computeCostAD start" << std::endl;
  std::shared_ptr<const Increment_> dx = std::dynamic_pointer_cast<const Increment_>(pv);
  ftlad_->setupAD(dx);
  ppad.enrollProcessor(ftlad_);
  Log::trace() << "CostJcDFI::computeCostAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
std::unique_ptr<GeneralizedDepartures>
CostJcDFI<MODEL, OBS>::multiplyCovar(const GeneralizedDepartures & dv1) const {
  const Increment_ & dx1 = dynamic_cast<const Increment_ &>(dv1);
  std::unique_ptr<Increment_> dx2(new Increment_(dx1));
  const double za = 1.0/alpha_;
  *dx2 *= za;
  return std::move(dx2);
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
std::unique_ptr<GeneralizedDepartures>
CostJcDFI<MODEL, OBS>::multiplyCoInv(const GeneralizedDepartures & dv1) const {
  const Increment_ & dx1 = dynamic_cast<const Increment_ &>(dv1);
  std::unique_ptr<Increment_> dx2(new Increment_(dx1));
  *dx2 *= alpha_;
  return std::move(dx2);
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
std::unique_ptr<GeneralizedDepartures> CostJcDFI<MODEL, OBS>::newDualVector() const {
  std::unique_ptr<Increment_> dx(new Increment_(*tlres_, vars_, vt_));
  return std::move(dx);
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
std::unique_ptr<GeneralizedDepartures> CostJcDFI<MODEL, OBS>::newGradientFG() const {
  return std::unique_ptr<Increment_>(new Increment_(*gradFG_));
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJcDFI<MODEL, OBS>::resetLinearization() {
  gradFG_.reset();
  ftlad_.reset();
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJCDFI_H_
