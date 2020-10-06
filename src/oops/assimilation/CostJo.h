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

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/Departures.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/Observers.h"
#include "oops/base/ObserversTLAD.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostBase.h"
#include "oops/base/PostBaseTLAD.h"
#include "oops/base/QCData.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Jo Cost Function
/*!
 * The CostJo class encapsulates the Jo term of the cost function.
 * The Observers to be called during the model integration is managed
 * inside the CostJo class.
 */

template<typename MODEL, typename OBS> class CostJo : public CostTermBase<MODEL, OBS>,
                                        private boost::noncopyable {
  typedef ControlVariable<MODEL, OBS>   CtrlVar_;
  typedef ControlIncrement<MODEL, OBS>  CtrlInc_;
  typedef Departures<OBS>               Departures_;
  typedef Observations<OBS>             Observations_;
  typedef Geometry<MODEL>               Geometry_;
  typedef State<MODEL>                  State_;
  typedef Increment<MODEL>              Increment_;
  typedef ObsErrors<OBS>                ObsErrors_;
  typedef ObsSpaces<OBS>                ObsSpaces_;
  typedef Observers<MODEL, OBS>         Observers_;
  typedef ObserversTLAD<MODEL, OBS>     ObserversTLAD_;
  typedef PostBaseTLAD<MODEL>           PostBaseTLAD_;
  typedef QCData<OBS>                   QCData_;

 public:
  /// Construct \f$ J_o\f$ from \f$ R\f$ and \f$ y_{obs}\f$.
  CostJo(const eckit::Configuration &, const eckit::mpi::Comm &,
         const util::DateTime &, const util::DateTime &,
         const eckit::mpi::Comm & ctime = oops::mpi::myself());

  /// Constructor added for generic 1d-var under development in ufo
  CostJo(const eckit::Configuration &, const util::DateTime &, const util::DateTime &,
         const ObsSpaces_ &);

  /// Destructor
  virtual ~CostJo() {}

  /// Initialize \f$ J_o\f$ before starting the integration of the model.
  std::shared_ptr<PostBase<State_> > initialize(const CtrlVar_ &,
                                                const eckit::Configuration &) override;
  std::shared_ptr<PostBaseTLAD_> initializeTraj(const CtrlVar_ &,
                                                const Geometry_ &,
                                                const eckit::Configuration &) override;
  /// Finalize \f$ J_o\f$ after the integration of the model.
  double finalize() override;
  void finalizeTraj() override;

  /// Initialize \f$ J_o\f$ before starting the TL run.
  std::shared_ptr<PostBaseTLAD_> setupTL(const CtrlInc_ &) const override;

  /// Initialize \f$ J_o\f$ before starting the AD run.
  std::shared_ptr<PostBaseTLAD_> setupAD(
           std::shared_ptr<const GeneralizedDepartures>, CtrlInc_ &) const override;

  /// Multiply by \f$ R\f$ and \f$ R^{-1}\f$.
  std::unique_ptr<GeneralizedDepartures>
    multiplyCovar(const GeneralizedDepartures &) const override;
  std::unique_ptr<GeneralizedDepartures>
    multiplyCoInv(const GeneralizedDepartures &) const override;

  /// Provide new departure.
  std::unique_ptr<GeneralizedDepartures> newDualVector() const override;

  /// Return gradient at first guess ie \f$ R^{-1} {\cal H}(x^t ) - y\f$.
  std::unique_ptr<GeneralizedDepartures> newGradientFG() const override;

  /// Reset obs operator trajectory.
  void resetLinearization() override;

  /// Print Jo
  double printJo(const Departures_ &, const Departures_ &) const;
  const ObsSpaces_ & obspaces() const {return obspace_;}

 private:
  eckit::LocalConfiguration obsconf_;
  ObsSpaces_ obspace_;
  Observations_ yobs_;
  std::unique_ptr<ObsErrors_> Rmat_;

  /// Configuration for current initialize/finalize pair
  std::unique_ptr<eckit::LocalConfiguration> currentConf_;

  /// Gradient at first guess : \f$ R^{-1} (H(x_{fg})-y_{obs}) \f$.
  std::unique_ptr<Departures_> gradFG_;

  /// Observers passed by \f$ J_o\f$ to the model during integration.
  mutable std::shared_ptr<Observers_> pobs_;

  /// Linearized observation operators.
  std::shared_ptr<ObserversTLAD_> pobstlad_;

  /// Storage for QC flags and obs error
  QCData_ qc_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostJo<MODEL, OBS>::CostJo(const eckit::Configuration & joConf, const eckit::mpi::Comm & comm,
                           const util::DateTime & winbgn, const util::DateTime & winend,
                           const eckit::mpi::Comm & ctime)
  : obsconf_(joConf), obspace_(obsconf_, comm, winbgn, winend, ctime),
    yobs_(obspace_, "ObsValue"),
    Rmat_(), currentConf_(), gradFG_(), pobs_(),
    pobstlad_(), qc_(obspace_)
{
  Log::trace() << "CostJo::CostJo done" << std::endl;
}

// =============================================================================
/// Constructor added for generic 1d-var under development in ufo
template<typename MODEL, typename OBS>
CostJo<MODEL, OBS>::CostJo(const eckit::Configuration & joConf,
                           const util::DateTime & winbgn, const util::DateTime & winend,
                           const ObsSpaces_ & localobs)
  : obsconf_(joConf), obspace_(localobs),
    yobs_(obspace_, "ObsValue"),
    Rmat_(), currentConf_(), gradFG_(), pobs_(),
    pobstlad_(), qc_(obspaces)
{
  Log::trace() << "CostJo::CostJo using local obs spaces done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
std::shared_ptr<PostBase<State<MODEL> > >
CostJo<MODEL, OBS>::initialize(const CtrlVar_ & xx, const eckit::Configuration & conf) {
  Log::trace() << "CostJo::initialize start" << std::endl;

  currentConf_.reset(new eckit::LocalConfiguration(conf));
  const int iterout = currentConf_->getInt("iteration");
  obsconf_.set("iteration", iterout);
  pobs_.reset(new Observers_(obsconf_, obspace_, xx.obsVar(), qc_));
  Log::trace() << "CostJo::initialize done" << std::endl;
  return pobs_;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
double CostJo<MODEL, OBS>::finalize() {
  Log::trace() << "CostJo::finalize start" << std::endl;
  const Observations_ & yeqv = pobs_->hofx();
  Log::info() << "Jo Observation Equivalent:" << std::endl << yeqv
              << "End Jo Observation Equivalent" << std::endl;
  const int iterout = currentConf_->getInt("iteration");

// Sace current QC flags and obs error
  const std::string obsname = "hofx" + std::to_string(iterout);
  const std::string qcname  = "EffectiveQC" + std::to_string(iterout);
  const std::string errname = "EffectiveError" + std::to_string(iterout);
  yeqv.save(obsname);
  for (size_t jj = 0; jj < obspace_.size(); ++jj) {
    qc_.obsErrors(jj)->mask(*qc_.qcFlags(jj));
    qc_.qcFlags(jj)->save(qcname);
    qc_.obsErrors(jj)->save(errname);
    qc_.obsErrors(jj)->save("EffectiveError");  // Obs error covariance is looking for that for now
  }

// Set observation error covariance
  Rmat_.reset(new ObsErrors_(obsconf_, obspace_));

// Perturb observations according to obs error statistics
  bool obspert = currentConf_->getBool("obs perturbations", false);
  if (obspert) {
    yobs_.perturb(*Rmat_);
    Log::info() << "Perturbed observations: " << yobs_ << std::endl;
  }

// Compute departures
  Departures_ ydep(yeqv - yobs_);
  Log::info() << "Jo Departures:" << std::endl << ydep << "End Jo Departures" << std::endl;

// Apply bias correction
  Departures_ bias(obspace_, "ObsBias", false);
  ydep += bias;
  Log::info() << "Jo Bias Corrected Departures:" << std::endl << ydep
          << "End Jo Bias Corrected Departures" << std::endl;

// Compute Jo
  if (!gradFG_) {
    gradFG_.reset(new Departures_(ydep));
  } else {
    *gradFG_ = ydep;
  }
  Rmat_->inverseMultiply(*gradFG_);

  double zjo = this->printJo(ydep, *gradFG_);

  if (currentConf_->has("diagnostics.departures")) {
    const std::string depname = currentConf_->getString("diagnostics.departures");
    ydep.save(depname);
  }

  pobs_.reset();
  currentConf_.reset();
  Log::trace() << "CostJo::finalize done" << std::endl;
  return zjo;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
std::shared_ptr<PostBaseTLAD<MODEL> >
CostJo<MODEL, OBS>::initializeTraj(const CtrlVar_ & xx, const Geometry_ &,
                              const eckit::Configuration & conf) {
  Log::trace() << "CostJo::initializeTraj start" << std::endl;
  pobstlad_.reset(new ObserversTLAD_(obsconf_, obspace_, xx.obsVar()));
  Log::trace() << "CostJo::initializeTraj done" << std::endl;
  return pobstlad_;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJo<MODEL, OBS>::finalizeTraj() {
  Log::trace() << "CostJo::finalizeTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
std::shared_ptr<PostBaseTLAD<MODEL> > CostJo<MODEL, OBS>::setupTL(const CtrlInc_ & dx) const {
  Log::trace() << "CostJo::setupTL start" << std::endl;
  ASSERT(pobstlad_);
  pobstlad_->setupTL(dx.obsVar());
  Log::trace() << "CostJo::setupTL done" << std::endl;
  return pobstlad_;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
std::shared_ptr<PostBaseTLAD<MODEL> > CostJo<MODEL, OBS>::setupAD(
                               std::shared_ptr<const GeneralizedDepartures> pv,
                               CtrlInc_ & dx) const {
  Log::trace() << "CostJo::setupAD start" << std::endl;
  ASSERT(pobstlad_);
  std::shared_ptr<const Departures_> dy = std::dynamic_pointer_cast<const Departures_>(pv);
  pobstlad_->setupAD(dy, dx.obsVar());
  Log::trace() << "CostJo::setupAD done" << std::endl;
  return pobstlad_;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
std::unique_ptr<GeneralizedDepartures>
CostJo<MODEL, OBS>::multiplyCovar(const GeneralizedDepartures & v1) const {
  Log::trace() << "CostJo::multiplyCovar start" << std::endl;
  std::unique_ptr<Departures_> y1(new Departures_(dynamic_cast<const Departures_ &>(v1)));
  Rmat_->multiply(*y1);
  return std::move(y1);
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
std::unique_ptr<GeneralizedDepartures>
CostJo<MODEL, OBS>::multiplyCoInv(const GeneralizedDepartures & v1) const {
  Log::trace() << "CostJo::multiplyCoInv start" << std::endl;
  std::unique_ptr<Departures_> y1(new Departures_(dynamic_cast<const Departures_ &>(v1)));
  Rmat_->inverseMultiply(*y1);
  return std::move(y1);
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
std::unique_ptr<GeneralizedDepartures> CostJo<MODEL, OBS>::newDualVector() const {
  Log::trace() << "CostJo::newDualVector start" << std::endl;
  std::unique_ptr<Departures_> ydep(new Departures_(obspace_));
  ydep->zero();
  Log::trace() << "CostJo::newDualVector done" << std::endl;
  return std::move(ydep);
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
std::unique_ptr<GeneralizedDepartures> CostJo<MODEL, OBS>::newGradientFG() const {
  return std::unique_ptr<Departures_>(new Departures_(*gradFG_));
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJo<MODEL, OBS>::resetLinearization() {
  Log::trace() << "CostJo::resetLinearization start" << std::endl;
  pobstlad_.reset();
  Log::trace() << "CostJo::resetLinearization done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
double CostJo<MODEL, OBS>::printJo(const Departures_ & dy, const Departures_ & grad) const {
  Log::trace() << "CostJo::printJo start" << std::endl;
  obspace_.printJo(dy, grad);

  double zjo = 0.0;
  for (std::size_t jj = 0; jj < dy.size(); ++jj) {
    const double zz = 0.5 * dot_product(dy[jj], grad[jj]);
    const unsigned nobs = grad[jj].nobs();
    if (nobs > 0) {
      Log::test() << "CostJo   : Nonlinear Jo(" << obspace_[jj].obsname() << ") = "
                  << zz << ", nobs = " << nobs << ", Jo/n = " << zz/nobs
                  << ", err = " << (*Rmat_)[jj].getRMSE() << std::endl;
    } else {
      Log::test() << "CostJo   : Nonlinear Jo(" << obspace_[jj].obsname() << ") = "
                  << zz << " --- No Observations" << std::endl;
      Log::warning() << "CostJo: No Observations!!!" << std::endl;
    }
    zjo += zz;
  }

  Log::trace() << "CostJo::printJo done" << std::endl;
  return zjo;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJO_H_
