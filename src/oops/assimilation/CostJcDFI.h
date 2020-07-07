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

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/DolphChebyshev.h"
#include "oops/base/PostBase.h"
#include "oops/base/PostBaseTLAD.h"
#include "oops/base/Variables.h"
#include "oops/base/WeightedDiff.h"
#include "oops/base/WeightedDiffTLAD.h"
#include "oops/base/WeightingFct.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

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
  typedef PostBaseTLAD<MODEL>           PostBaseTLAD_;
  typedef State<MODEL>                  State_;

 public:
/// Construct \f$ J_c\f$.
  CostJcDFI(const eckit::Configuration &, const Geometry_ &, const util::DateTime &,
            const util::Duration &, const util::Duration & tstep = util::Duration(0));

/// Destructor
  virtual ~CostJcDFI() {}

/// Initialize before nonlinear model integration.
  std::shared_ptr<PostBase<State_> > initialize(const CtrlVar_ &,
                                                const eckit::Configuration &) override;
  std::shared_ptr<PostBaseTLAD<MODEL> > initializeTraj(const CtrlVar_ &, const Geometry_ &,
                                                       const eckit::Configuration &) override;

/// Finalize computation after nonlinear model integration.
  double finalize() override;
  void finalizeTraj() override;

/// Initialize \f$ J_c\f$ before starting the TL run.
  std::shared_ptr<PostBaseTLAD_> setupTL(const CtrlInc_ &) const override;

/// Initialize \f$ J_c\f$ before starting the AD run.
  std::shared_ptr<PostBaseTLAD_> setupAD(
           std::shared_ptr<const GeneralizedDepartures>, CtrlInc_ &) const override;

/// Multiply by \f$ C\f$ and \f$ C^{-1}\f$.
  Increment_ * multiplyCovar(const GeneralizedDepartures &) const override;
  Increment_ * multiplyCoInv(const GeneralizedDepartures &) const override;

/// Provide new increment.
  Increment_ * newDualVector() const override;

/// Gradient of \f$ J_c\f$ at first guess.
  Increment_ * newGradientFG() const override {return new Increment_(*gradFG_);}

/// Reset trajectory.
  void resetLinearization() override;

 private:
  const eckit::LocalConfiguration conf_;
  util::DateTime vt_;
  util::Duration span_;
  double alpha_;
  std::unique_ptr<WeightingFct> wfct_;
  std::unique_ptr<Increment_> gradFG_;
  const Geometry_ resol_;
  const util::Duration tstep_;
  std::unique_ptr<Geometry_> tlres_;
  util::Duration tlstep_;
  mutable std::shared_ptr<WeightedDiff<MODEL, Increment_, State_> > filter_;
  mutable std::shared_ptr<WeightedDiffTLAD<MODEL> > ftlad_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostJcDFI<MODEL, OBS>::CostJcDFI(const eckit::Configuration & conf, const Geometry_ & resol,
                            const util::DateTime & vt, const util::Duration & span,
                            const util::Duration & tstep)
  : conf_(conf), vt_(vt), span_(span), alpha_(0), wfct_(), gradFG_(),
    resol_(resol), tstep_(tstep), tlres_(), tlstep_(), filter_()
{
  alpha_ = conf.getDouble("alpha");
  if (conf.has("ftime")) vt_ = util::DateTime(conf.getString("ftime"));
  if (conf.has("span")) span_ = util::Duration(conf.getString("span"));
//  wfct_.reset(WeightFactory::create(config)); YT
  wfct_.reset(new DolphChebyshev(conf_));
  Log::debug() << "CostJcDFI created vt = " << vt_ << ", span = " << span_ << std::endl;
  Log::trace() << "CostJcDFI created" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
std::shared_ptr<PostBase<State<MODEL> > >
CostJcDFI<MODEL, OBS>::initialize(const CtrlVar_ &, const eckit::Configuration &) {
  filter_.reset(new WeightedDiff<MODEL, Increment_, State_>(conf_, vt_, span_,
                                                            tstep_, resol_, *wfct_));
  return filter_;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
double CostJcDFI<MODEL, OBS>::finalize() {
  double zz = 0.5 * alpha_;
  std::unique_ptr<Increment_> dx(filter_->releaseDiff());
  zz *= dot_product(*dx, *dx);
  Log::test() << "CostJcDFI: Nonlinear Jc = " << zz << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
std::shared_ptr<PostBaseTLAD<MODEL> >
CostJcDFI<MODEL, OBS>::initializeTraj(const CtrlVar_ &, const Geometry_ & tlres,
                                 const eckit::Configuration & innerConf) {
  tlres_.reset(new Geometry_(tlres));
  tlstep_ = util::Duration(innerConf.getString("linearmodel.tstep"));
  ftlad_.reset(new WeightedDiffTLAD<MODEL>(conf_, vt_, span_, tstep_, *tlres_, *wfct_));
  return ftlad_;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJcDFI<MODEL, OBS>::finalizeTraj() {
  gradFG_.reset(ftlad_->releaseDiff());
  *gradFG_ *= alpha_;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
Increment<MODEL> * CostJcDFI<MODEL, OBS>::newDualVector() const {
  const Variables vars(conf_);
  Increment_ * dx = new Increment_(*tlres_, vars, vt_);
  return dx;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
std::shared_ptr<PostBaseTLAD<MODEL> >
CostJcDFI<MODEL, OBS>::setupTL(const CtrlInc_ &) const {
  ftlad_->setupTL(*tlres_);
  return ftlad_;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
std::shared_ptr<PostBaseTLAD<MODEL> >
CostJcDFI<MODEL, OBS>::setupAD(std::shared_ptr<const GeneralizedDepartures> pv,
                          CtrlInc_ &) const {
  std::shared_ptr<const Increment_>
    dx = std::dynamic_pointer_cast<const Increment_>(pv);
  ftlad_->setupAD(dx);
  return ftlad_;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
Increment<MODEL> * CostJcDFI<MODEL, OBS>::multiplyCovar(const GeneralizedDepartures & dv1) const {
  const Increment_ & dx1 = dynamic_cast<const Increment_ &>(dv1);
  Increment_ * dx2 = new Increment_(dx1);
  const double za = 1.0/alpha_;
  *dx2 *= za;
  return dx2;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
Increment<MODEL> * CostJcDFI<MODEL, OBS>::multiplyCoInv(const GeneralizedDepartures & dv1) const {
  const Increment_ & dx1 = dynamic_cast<const Increment_ &>(dv1);
  Increment_ * dx2 = new Increment_(dx1);
  *dx2 *= alpha_;
  return dx2;
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
