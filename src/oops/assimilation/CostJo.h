/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_COSTJO_H_
#define OOPS_ASSIMILATION_COSTJO_H_

#include <string>
#include <boost/noncopyable.hpp>
#include <boost/pointer_cast.hpp>
#include <boost/scoped_ptr.hpp>
#include <boost/shared_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/Departures.h"
#include "oops/base/LinearObsOperators.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/Observer.h"
#include "oops/base/ObserverTLAD.h"
#include "oops/base/ObsFilters.h"
#include "oops/base/ObsOperators.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostBase.h"
#include "oops/base/PostBaseTLAD.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/ObsAuxIncrement.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Jo Cost Function
/*!
 * The CostJo class encapsulates the Jo term of the cost function.
 * The Observer to be called during the model integration is managed
 * inside the CostJo class.
 */

template<typename MODEL> class CostJo : public CostTermBase<MODEL>,
                                        private boost::noncopyable {
  typedef ControlVariable<MODEL>     CtrlVar_;
  typedef ControlIncrement<MODEL>    CtrlInc_;
  typedef Departures<MODEL>          Departures_;
  typedef Observations<MODEL>        Observations_;
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;
  typedef Increment<MODEL>           Increment_;
  typedef ObsAuxIncrement<MODEL>     ObsAuxIncr_;
  typedef ObsFilters<MODEL>          ObsFilters_;
  typedef ObsOperators<MODEL>        ObsOperator_;
  typedef ObsSpaces<MODEL>           ObsSpace_;
  typedef PostBaseTLAD<MODEL>        PostBaseTLAD_;
  typedef LinearObsOperators<MODEL>  LinearObsOperator_;

 public:
  /// Construct \f$ J_o\f$ from \f$ R\f$ and \f$ y_{obs}\f$.
  CostJo(const eckit::Configuration &, const util::DateTime &, const util::DateTime &,
         const util::Duration &, const bool subwindows = false);

  /// Destructor
  virtual ~CostJo() {}

  /// Initialize \f$ J_o\f$ before starting the integration of the model.
  boost::shared_ptr<PostBase<State_> > initialize(const CtrlVar_ &) const override;
  boost::shared_ptr<PostBase<State_> > initializeTraj(const CtrlVar_ &,
                                                      const Geometry_ &,
                                                      const eckit::Configuration &) override;
  /// Finalize \f$ J_o\f$ after the integration of the model.
  double finalize(const eckit::Configuration &) const override;
  double finalizeTraj(const eckit::Configuration &) override;

  /// Initialize \f$ J_o\f$ before starting the TL run.
  boost::shared_ptr<PostBaseTLAD_> setupTL(const CtrlInc_ &) const override;

  /// Initialize \f$ J_o\f$ before starting the AD run.
  boost::shared_ptr<PostBaseTLAD_> setupAD(
           boost::shared_ptr<const GeneralizedDepartures>, CtrlInc_ &) const override;

  /// Multiply by \f$ R\f$ and \f$ R^{-1}\f$.
  Departures_ * multiplyCovar(const GeneralizedDepartures &) const override;
  Departures_ * multiplyCoInv(const GeneralizedDepartures &) const override;

  /// Provide new departure.
  Departures_ * newDualVector() const override;

  /// Return gradient at first guess ie \f$ R^{-1} {\cal H}(x^t ) - y\f$.
  Departures_ * newGradientFG() const override {return new Departures_(*gradFG_);}

  /// Reset obs operator trajectory.
  void resetLinearization() override;

  /// Print Jo
  double printJo(const Departures_ &, const Departures_ &) const;

 private:
  const ObsSpace_ obspace_;
  const ObsOperator_ hop_;
  Observations_ yobs_;
  ObsErrors<MODEL> R_;

  /// Gradient at first guess : \f$ R^{-1} (H(x_{fg})-y_{obs}) \f$.
  boost::scoped_ptr<Departures_> gradFG_;

  /// Observer passed by \f$ J_o\f$ to the model during integration.
  mutable boost::shared_ptr<Observer<MODEL, State_> > pobs_;

  /// Time slot.
  const util::Duration tslot_;

  /// Linearized observation operators.
  boost::shared_ptr<LinearObsOperator_> hoptlad_;
  const bool subwindows_;
  bool ltraj_;
};

// =============================================================================

template<typename MODEL>
CostJo<MODEL>::CostJo(const eckit::Configuration & joConf,
                      const util::DateTime & winbgn, const util::DateTime & winend,
                      const util::Duration & tslot, const bool subwindows)
  : obspace_(joConf, winbgn, winend),
    hop_(obspace_), yobs_(obspace_), R_(obspace_),
    gradFG_(), pobs_(), tslot_(tslot),
    hoptlad_(), subwindows_(subwindows), ltraj_(false)
{
  Log::debug() << "CostJo:setup tslot_ = " << tslot_ << std::endl;
  yobs_.read(joConf);
  Log::trace() << "CostJo:CostJo done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
boost::shared_ptr<PostBase<State<MODEL> > >
CostJo<MODEL>::initialize(const CtrlVar_ & xx) const {
  ASSERT(ltraj_ == false);
  const eckit::LocalConfiguration conf;
  ObsFilters_ filter(obspace_, conf);
  pobs_.reset(new Observer<MODEL, State_>(obspace_, hop_, xx.obsVar(), filter,
                                          tslot_, subwindows_));
  return pobs_;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double CostJo<MODEL>::finalize(const eckit::Configuration & conf) const {
  ASSERT(ltraj_ == false);

  boost::scoped_ptr<Observations_> yeqv(pobs_->release());
  Log::info() << "Jo Observation Equivalent:" << *yeqv << std::endl;

  Departures_ ydep(*yeqv - yobs_);
  Log::info() << "Jo Departures:" << ydep << std::endl;

  boost::scoped_ptr<Departures_> grad(R_.inverseMultiply(ydep));

  double zjo = this->printJo(ydep, *grad);

  if (conf.has("departures")) {
    const std::string depname = conf.getString("departures");
    ydep.save(depname);
  }

  pobs_.reset();
  return zjo;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
boost::shared_ptr<PostBase<State<MODEL> > >
CostJo<MODEL>::initializeTraj(const CtrlVar_ & xx, const Geometry_ &,
                              const eckit::Configuration & conf) {
  ltraj_ = true;
  hoptlad_.reset(new LinearObsOperator_(obspace_));
  ObsFilters_ filter(obspace_, conf);
  pobs_.reset(new Observer<MODEL, State_>(obspace_, hop_, xx.obsVar(), filter,
                                          tslot_, subwindows_, hoptlad_));
  return pobs_;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double CostJo<MODEL>::finalizeTraj(const eckit::Configuration & conf) {
  ASSERT(ltraj_);

  boost::scoped_ptr<Observations_> yeqv(pobs_->release());
  Log::info() << "Jo Observation Equivalent:" << *yeqv << std::endl;
  R_.linearize(*yeqv);

  Departures_ ydep(*yeqv - yobs_);
  Log::info() << "Jo Departures:" << ydep << std::endl;

  gradFG_.reset(R_.inverseMultiply(ydep));

  double zjo = this->printJo(ydep, *gradFG_);

  if (conf.has("departures")) {
    const std::string depname = conf.getString("departures");
    ydep.save(depname);
  }

  ltraj_ = false;

  pobs_.reset();
  return zjo;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
boost::shared_ptr<PostBaseTLAD<MODEL> > CostJo<MODEL>::setupTL(const CtrlInc_ & dx) const {
  ASSERT(hoptlad_);
  boost::shared_ptr<PostBaseTLAD_> spobs;
  spobs.reset(new ObserverTLAD<MODEL>(obspace_, *hoptlad_, dx.obsVar(),
                                      tslot_, subwindows_));
  return spobs;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
boost::shared_ptr<PostBaseTLAD<MODEL> > CostJo<MODEL>::setupAD(
                               boost::shared_ptr<const GeneralizedDepartures> pv,
                               CtrlInc_ & dx) const {
  ASSERT(hoptlad_);
  boost::shared_ptr<const Departures_> dy = boost::dynamic_pointer_cast<const Departures_>(pv);
  boost::shared_ptr<PostBaseTLAD_> spobs;
  spobs.reset(new ObserverTLAD<MODEL>(obspace_, *hoptlad_, dy, dx.obsVar(),
                                      tslot_, subwindows_));
  return spobs;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Departures<MODEL> * CostJo<MODEL>::multiplyCovar(const GeneralizedDepartures & v1) const {
  const Departures_ & y1 = dynamic_cast<const Departures_ &>(v1);
  return R_.multiply(y1);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Departures<MODEL> * CostJo<MODEL>::multiplyCoInv(const GeneralizedDepartures & v1) const {
  const Departures_ & y1 = dynamic_cast<const Departures_ &>(v1);
  return R_.inverseMultiply(y1);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Departures<MODEL> * CostJo<MODEL>::newDualVector() const {
  Departures_ * ydep = new Departures_(obspace_);
  ydep->zero();
  return ydep;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void CostJo<MODEL>::resetLinearization() {
  hoptlad_.reset();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double CostJo<MODEL>::printJo(const Departures_ & dy, const Departures_ & grad) const {
  obspace_.printJo(dy, grad);

  double zjo = 0.0;
  for (std::size_t jj = 0; jj < dy.size(); ++jj) {
    const double zz = 0.5 * dot_product(dy[jj], grad[jj]);
    const unsigned nobs = dy[jj].size();
    if (nobs > 0) {
      Log::test() << "CostJo   : Nonlinear Jo = " << zz
                  << ", nobs = " << nobs << ", Jo/n = " << zz/nobs
                  << ", err = " << R_[jj].getRMSE() << std::endl;
    } else {
      Log::test() << "CostJo   : Nonlinear Jo = " << zz << " --- No Observations" << std::endl;
      Log::warning() << "CostJo: No Observations!!!" << std::endl;
    }
    zjo += zz;
  }

  return zjo;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJO_H_
