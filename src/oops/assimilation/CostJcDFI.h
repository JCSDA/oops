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

#include <boost/pointer_cast.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "util/Logger.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/DolphChebyshev.h"
#include "oops/base/PostBase.h"
#include "oops/base/PostBaseTL.h"
#include "oops/base/PostBaseAD.h"
#include "oops/base/WeightedDiff.h"
#include "oops/base/WeightedDiffTL.h"
#include "oops/base/WeightedDiffAD.h"
#include "oops/base/WeightingFct.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/interface/Variables.h"
#include "util/DateTime.h"
#include "util/Duration.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Jc DFI Cost Function
/*!
 * Digital filter based constraint term for the cost function.
 */

template<typename MODEL> class CostJcDFI : public CostTermBase<MODEL> {
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef ControlVariable<MODEL>     CtrlVar_;
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;
  typedef Increment<MODEL>           Increment_;
  typedef Variables<MODEL>           Variables_;

 public:
/// Construct \f$ J_c\f$.
  CostJcDFI(const eckit::Configuration &, const Geometry_ &, const util::DateTime &,
            const util::Duration &, const util::Duration & tstep = util::Duration(0));

/// Destructor
  virtual ~CostJcDFI() {}

/// Initialize before nonlinear model integration.
  boost::shared_ptr<PostBase<State_> > initialize(const CtrlVar_ &) const override;
  boost::shared_ptr<PostBase<State_> > initializeTraj(const CtrlVar_ &,
                                                      const Geometry_ &,
                                                      const eckit::Configuration &) override;

/// Finalize computation after nonlinear model integration.
  double finalize(const eckit::Configuration &) const override;
  double finalizeTraj(const eckit::Configuration &) override;

/// Initialize \f$ J_c\f$ before starting the TL run.
  boost::shared_ptr<PostBaseTL<Increment_> > setupTL(const CtrlInc_ &) const override;

/// Initialize \f$ J_c\f$ before starting the AD run.
  boost::shared_ptr<PostBaseAD<Increment_> > setupAD(
           boost::shared_ptr<const GeneralizedDepartures>, CtrlInc_ &) const override;

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
  boost::scoped_ptr<WeightingFct> wfct_;
  boost::scoped_ptr<Increment_> gradFG_;
  const Geometry_ resol_;
  const util::Duration tstep_;
  boost::scoped_ptr<Geometry_> tlres_;
  util::Duration tlstep_;
  mutable boost::shared_ptr<WeightedDiff<MODEL, Increment_, State_> > filter_;
  bool ltraj_;
};

// =============================================================================

template<typename MODEL>
CostJcDFI<MODEL>::CostJcDFI(const eckit::Configuration & conf, const Geometry_ & resol,
                            const util::DateTime & vt, const util::Duration & span,
                            const util::Duration & tstep)
  : conf_(conf), vt_(vt), span_(span), alpha_(0), wfct_(), gradFG_(),
    resol_(resol), tstep_(tstep), tlres_(), tlstep_(), filter_(), ltraj_(false)
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

template<typename MODEL>
boost::shared_ptr<PostBase<State<MODEL> > >
CostJcDFI<MODEL>::initialize(const CtrlVar_ &) const {
  ASSERT(ltraj_ == false);
  filter_.reset(new WeightedDiff<MODEL, Increment_, State_>(vt_, span_, resol_,
                                                            conf_, tstep_, *wfct_));
  return filter_;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double CostJcDFI<MODEL>::finalize(const eckit::Configuration &) const {
  ASSERT(ltraj_ == false);
  double zz = 0.5 * alpha_;
  boost::scoped_ptr<Increment_> dx(filter_->releaseDiff());
  zz *= dot_product(*dx, *dx);
  Log::test() << "CostJcDFI: Nonlinear Jc = " << zz << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
boost::shared_ptr<PostBase<State<MODEL> > >
CostJcDFI<MODEL>::initializeTraj(const CtrlVar_ &, const Geometry_ & tlres,
                                 const eckit::Configuration & innerConf) {
  ltraj_ = true;
  tlres_.reset(new Geometry_(tlres));
  tlstep_ = util::Duration(innerConf.getString("linearmodel.tstep"));
  filter_.reset(new WeightedDiff<MODEL, Increment_, State_>(vt_, span_, resol_,
                                                            conf_, tstep_, *wfct_));
  return filter_;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double CostJcDFI<MODEL>::finalizeTraj(const eckit::Configuration &) {
  ASSERT(ltraj_ == true);
  double zz = 0.5 * alpha_;
  gradFG_.reset(filter_->releaseDiff());
  zz *= dot_product(*gradFG_, *gradFG_);
  *gradFG_ *= alpha_;
  ltraj_ = false;
  Log::test() << "CostJcDFI: Nonlinear Jc = " << zz << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> * CostJcDFI<MODEL>::newDualVector() const {
  Variables_ vars(conf_);
  Increment_ * dx = new Increment_(*tlres_, vars, vt_);
  return dx;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
boost::shared_ptr<PostBaseTL<Increment<MODEL> > >
CostJcDFI<MODEL>::setupTL(const CtrlInc_ &) const {
  boost::shared_ptr<WeightedDiffTL<MODEL, Increment_> > filterTL(
    new WeightedDiffTL<MODEL, Increment_>(vt_, span_, *tlres_, conf_, tlstep_, *wfct_));
  return filterTL;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
boost::shared_ptr<PostBaseAD<Increment<MODEL> > >
CostJcDFI<MODEL>::setupAD(boost::shared_ptr<const GeneralizedDepartures> pv,
                          CtrlInc_ &) const {
  boost::shared_ptr<const Increment_>
    dx = boost::dynamic_pointer_cast<const Increment_>(pv);
  boost::shared_ptr<WeightedDiffAD<Increment_> > filterAD(
    new WeightedDiffAD<Increment_>(vt_, span_, tlstep_, *wfct_, dx));
  return filterAD;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> * CostJcDFI<MODEL>::multiplyCovar(const GeneralizedDepartures & dv1) const {
  const Increment_ & dx1 = dynamic_cast<const Increment_ &>(dv1);
  Increment_ * dx2 = new Increment_(dx1);
  const double za = 1.0/alpha_;
  *dx2 *= za;
  return dx2;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Increment<MODEL> * CostJcDFI<MODEL>::multiplyCoInv(const GeneralizedDepartures & dv1) const {
  const Increment_ & dx1 = dynamic_cast<const Increment_ &>(dv1);
  Increment_ * dx2 = new Increment_(dx1);
  *dx2 *= alpha_;
  return dx2;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJcDFI<MODEL>::resetLinearization() {
  gradFG_.reset();
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJCDFI_H_
