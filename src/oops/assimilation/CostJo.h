/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Crown Copyright 2023, the Met Office.
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
#include "eckit/exception/Exceptions.h"
#include "oops/assimilation/ControlIncrement.h"
#include "oops/assimilation/ControlVariable.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/Departures.h"
#include "oops/base/Geometry.h"
#include "oops/base/GetValueTLADs.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/Observers.h"
#include "oops/base/ObserversTLAD.h"
#include "oops/base/ObsOperatorBase.h"
#include "oops/base/ObsOperatorPert.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/base/State.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

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
  typedef GetValueTLADs<MODEL, OBS>     GetValueTLADs_;
  typedef State<MODEL>                  State_;
  typedef ObsDataVector<OBS, int>       ObsDataInt_;
  typedef ObsErrors<OBS>                ObsErrors_;
  typedef ObsSpaces<OBS>                ObsSpaces_;
  typedef Observers<MODEL, OBS>         Observers_;
  typedef ObserversTLAD<MODEL, OBS>     ObserversTLAD_;
  typedef ObsOperatorBase<OBS>          ObsOperatorBase_;
  typedef ObsOperatorPert<OBS>          ObsOperatorPert_;
  typedef PostProcessor<State_>         PostProc_;
  typedef PostProcessorTLAD<MODEL>      PostProcTLAD_;

 public:
  /// Construct \f$ J_o\f$ from \f$ R\f$ and \f$ y_{obs}\f$.
  CostJo(const eckit::Configuration &, const eckit::mpi::Comm &,
         const util::TimeWindow &,
         const eckit::mpi::Comm & ctime = oops::mpi::myself());

  /// Destructor
  virtual ~CostJo();

  /// Initialize \f$ J_o\f$ before starting the integration of the model.
  void setPostProc(const CtrlVar_ &, const eckit::Configuration &, PostProc_ &) override;
  /// Finalize \f$ J_o\f$ after the integration of the model.
  double computeCost() override;
  void printCostTestHack() override;

  /// Initialize \f$ J_o\f$ for the trajectory run
  void setPostProcTraj(const CtrlVar_ &, const eckit::Configuration &,
                       const Geometry_ &, PostProcTLAD_ &) override;
  void computeCostTraj() override;

  /// Initialize \f$ J_o\f$ before starting the TL run.
  void setPostProcTL(const CtrlInc_ &, PostProcTLAD_ &) const override;
  void computeCostTL(const CtrlInc_ &, GeneralizedDepartures &) const override;

  /// Initialize \f$ J_o\f$ before starting the AD run.
  void computeCostAD(std::shared_ptr<const GeneralizedDepartures>,
                     CtrlInc_ &, PostProcTLAD_ &) const override;
  void setPostProcAD() const override;

  /// Multiply by \f$ R\f$ and \f$ R^{-1}\f$.
  std::unique_ptr<GeneralizedDepartures>
    multiplyCovar(const GeneralizedDepartures &) const override;
  std::unique_ptr<GeneralizedDepartures>
    multiplyCoInv(const GeneralizedDepartures &) const override;

  /// Provide new departure.
  std::unique_ptr<GeneralizedDepartures> newDualVector() const override;

  /// Return gradient at first guess ie \f$ R^{-1} {\cal H}(x^t ) - y\f$.
  std::unique_ptr<GeneralizedDepartures> newGradientFG() const override;

  void setObsPert(const Geometry_ &, const Variables &);

  /// Reset obs operator trajectory.
  void resetLinearization() override;

  /// Append new observations
  void appendObs(const eckit::Configuration & appendConfig);

  /// Accessors
  const ObsSpaces_ & obspaces() const {return obspaces_;}

 private:
  double printJo(size_t, Departures_ &, std::ostream &) const;

  const eckit::LocalConfiguration conf_;
  ObsSpaces_ obspaces_;
  std::unique_ptr<Observations_> yobs_;
  std::unique_ptr<ObsErrors_> Rmat_;
  std::unique_ptr<Observers_> observers_;
  std::vector<ObsDataInt_> qcflags_;

  /// Jo Gradient at first guess : \f$ R^{-1} (H(x_{fg})-y_{obs}) \f$.
  std::unique_ptr<Departures_> gradFG_;

  /// Linearized observation operators.
  std::shared_ptr<ObserversTLAD_> obstlad_;

  /// Configuration for current initialize/finalize pair
  std::unique_ptr<eckit::LocalConfiguration> currentConf_;
  std::unique_ptr<Departures_> yhack_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostJo<MODEL, OBS>::CostJo(const eckit::Configuration & joConf, const eckit::mpi::Comm & comm,
                           const util::TimeWindow & timeWindow,
                           const eckit::mpi::Comm & ctime)
  : conf_(joConf),
    obspaces_(eckit::LocalConfiguration(joConf, "observers"), comm, timeWindow, ctime),
    yobs_(), Rmat_(), observers_(), gradFG_(), obstlad_(), currentConf_()
{
  observers_.reset(new Observers_(obspaces_, joConf));
  Rmat_.reset(new ObsErrors_(eckit::LocalConfiguration(joConf, "observers"), obspaces_));
  for (size_t jj = 0; jj < obspaces_.size(); ++jj) {
    ObsDataInt_ qcflags(obspaces_[jj], obspaces_[jj].obsvariables());
    qcflags_.push_back(qcflags);
  }
  Log::trace() << "CostJo::CostJo" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
CostJo<MODEL, OBS>::~CostJo() {
  obspaces_.save();
  Log::trace() << "CostJo::~CostJo" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJo<MODEL, OBS>::setPostProc(const CtrlVar_ & xx, const eckit::Configuration & conf,
                                     PostProc_ & pp) {
  Log::trace() << "CostJo::setPostProc start" << std::endl;
  gradFG_.reset();
  currentConf_.reset(new eckit::LocalConfiguration(conf));
  observers_->initialize(xx.state().geometry(), xx.obsVar(), *Rmat_, pp, conf);
  Log::trace() << "CostJo::setPostProc done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
double CostJo<MODEL, OBS>::computeCost() {
  Log::trace() << "CostJo::computeCost start" << std::endl;
  bool firstOuterLoop = !yobs_;
  bool lastOuterLoop =
    (currentConf_->getInt("iteration") == currentConf_->getInt("total iterations"));

  // Obs, simulated obs and departures (held here for nice prints and diagnostics)
  Observations_ yeqv(obspaces_);
  observers_->finalize(yeqv, qcflags_);

  if (firstOuterLoop) {
    yobs_.reset(new Observations_(obspaces_, "ObsValue"));
    Log::info() << "CostJo Observations: " << *yobs_ << std::endl;
    if (conf_.getBool("obs perturbations", false)) {
      // Perturb observations according to obs error statistics and save to output file
      yobs_->perturb(*Rmat_);
      Log::info() << "Perturbed observations: " << *yobs_ << std::endl;
      yobs_->save("EffectiveObsValue");
    }
  }

  // Compute observations departures and save to output file
  Departures_ ydep(yeqv - *yobs_);
  const std::string depname = firstOuterLoop ? "ombg"
                            : lastOuterLoop  ? "oman"
                            : "innov" + currentConf_->getString("iteration");
  Departures_ minusydep = ydep;
  minusydep *= (-1.0);
  minusydep.save(depname);

  // Gradient at first guess (to define inner loop rhs)
  gradFG_.reset(new Departures_(ydep));
  Rmat_->inverseMultiply(*gradFG_);

  // Print diagnostics
  Log::info() << "Jo Observations:" << std::endl << *yobs_
          << "End Jo Observations"  << std::endl;

  Log::info() << "Jo Observations Equivalent:" << std::endl << yeqv
          << "End Jo Observations Equivalent"  << std::endl;

  Log::info() << "Jo Bias Corrected Departures:" << std::endl << ydep
          << "End Jo Bias Corrected Departures"  << std::endl;

  Log::info() << "Jo Observations Errors:" << std::endl << *Rmat_
          << "End Jo Observations Errors"  << std::endl;

  // Print Jo table
  double zjo = 0.0;
  for (size_t jj = 0; jj < obspaces_.size(); ++jj) {
    zjo += this->printJo(jj, ydep, Log::info());
  }
// TEMPORARY HACK START
//  for (size_t jj = 0; jj < obspaces_.size(); ++jj) {
//    this->printJo(jj, ydep, Log::test());
//  }
  yhack_.reset(new Departures_(ydep));
// TEMPORARY HACK END

  Log::info() << "CostJo   : Nonlinear Jo = " << zjo << std::endl;

  Log::trace() << "CostJo::computeCost done" << std::endl;
  return zjo;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJo<MODEL, OBS>::printCostTestHack() {
  for (size_t jj = 0; jj < obspaces_.size(); ++jj) {
    this->printJo(jj, *yhack_, Log::test());
  }
  yhack_.reset();
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
double CostJo<MODEL, OBS>::printJo(size_t jj, Departures_ & ydep, std::ostream & os) const {
  Log::trace() << "CostJo::printJo start" << std::endl;

  double zz = 0.0;
  const unsigned nobs = (*gradFG_)[jj].nobs();
  os << "CostJo   : Nonlinear Jo(" << obspaces_[jj].obsname() << ") = ";

  if (nobs > 0) {
    zz = 0.5 * dot_product(ydep[jj], (*gradFG_)[jj]);
    const double err = (*Rmat_)[jj].getRMSE();
    os << zz << ", nobs = " << nobs << ", Jo/n = " << zz/nobs << ", err = " << err;
  } else {
    os << zz << " --- No Observations";
  }

  std::vector<eckit::LocalConfiguration> subconfs =
      eckit::LocalConfiguration(conf_, "observers").getSubConfigurations();
  Log::trace() << "CostJo::printJo subconf" << subconfs[jj] << std::endl;
  if (subconfs[jj].getBool("monitoring only", false)) {
    os << " (Monitoring only)";
    zz = 0.0;
  }
  os << std::endl;

  Log::trace() << "CostJo::printJo done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJo<MODEL, OBS>::setPostProcTraj(const CtrlVar_ & xx, const eckit::Configuration & conf,
                                         const Geometry_ & lowres, PostProcTLAD_ & pptraj) {
  Log::trace() << "CostJo::setPostProcTraj start" << std::endl;
  obstlad_.reset(new ObserversTLAD_(obspaces_, conf_));
  obstlad_->initializeTraj(lowres, xx.obsVar(), pptraj);
  Log::trace() << "CostJo::setPostProcTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJo<MODEL, OBS>::computeCostTraj() {
  obstlad_->finalizeTraj(qcflags_);
  Log::trace() << "CostJo::computeCostTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJo<MODEL, OBS>::setPostProcTL(const CtrlInc_ & dx, PostProcTLAD_ & pptl) const {
  Log::trace() << "CostJo::setPostProcTL start" << std::endl;
  ASSERT(obstlad_);
  obstlad_->initializeTL(pptl);
  Log::trace() << "CostJo::setPostProcTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJo<MODEL, OBS>::computeCostTL(const CtrlInc_ & dx, GeneralizedDepartures & gdep) const {
  Log::trace() << "CostJo::computeCostTL start" << std::endl;
  ASSERT(obstlad_);
  Departures_ & ydep = dynamic_cast<Departures_ &>(gdep);
  obstlad_->finalizeTL(dx.obsVar(), ydep);
  Log::trace() << "CostJo::computeCostTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJo<MODEL, OBS>::computeCostAD(std::shared_ptr<const GeneralizedDepartures> pv,
                                       CtrlInc_ & dx, PostProcTLAD_ & ppad) const {
  Log::trace() << "CostJo::computeCostAD start" << std::endl;
  ASSERT(obstlad_);
  std::shared_ptr<const Departures_> dy = std::dynamic_pointer_cast<const Departures_>(pv);
  obstlad_->initializeAD(*dy, dx.obsVar(), ppad);
  Log::trace() << "CostJo::computeCostAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJo<MODEL, OBS>::setPostProcAD() const {
  Log::trace() << "CostJo::setPostProcAD start" << std::endl;
  ASSERT(obstlad_);
  obstlad_->finalizeAD();
  Log::trace() << "CostJo::setPostProcAD done" << std::endl;
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
  std::unique_ptr<Departures_> ydep(new Departures_(obspaces_));
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
void CostJo<MODEL, OBS>::setObsPert(const Geometry_ & geom, const Variables & incVars) {
  std::vector<eckit::LocalConfiguration> subconfs =
      eckit::LocalConfiguration(conf_, "observers").getSubConfigurations();

//  Reset the observation operators to be the linear observation operators contained in *linObsTLAD,
//  wrapped inside ObsOperatorPert
  std::vector<std::unique_ptr<ObsOperatorBase_>> obsOpBases;
  for (std::size_t jj = 0; jj < obspaces_.size(); ++jj) {
    const eckit::LocalConfiguration obsconf(subconfs[jj], "obs operator");
    obsOpBases.push_back(std::make_unique<ObsOperatorPert_>(obspaces_[jj], obsconf,
                                                             (*obstlad_)[jj].linObsOp(), true));
  }
  observers_->resetObsPert(geom, std::move(obsOpBases), obstlad_->posts(), incVars);
  yobs_.reset(new Observations_(this->obspaces(), ""));
  // Perturb observations according to obs error statistics and save to output file
  yobs_->perturb(*Rmat_);
  Log::info() << "Perturbed observations: " << *yobs_ << std::endl;
  yobs_->save("EffectiveObsValue");
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJo<MODEL, OBS>::resetLinearization() {
  Log::trace() << "CostJo::resetLinearization start" << std::endl;
  obstlad_.reset();
  Log::trace() << "CostJo::resetLinearization done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJo<MODEL, OBS>::appendObs(const eckit::Configuration & appendConfig) {
  Log::trace() << "CostJo::appendObs start" << std::endl;

  if (conf_.getBool("obs perturbations", false)) {
    throw eckit::NotImplemented("CostJo: appending to perturbed obs not implemented yet", Here());
  }

  if (!yobs_) {
    ABORT("CostJo::appendObs: Should not append obs in the first outer loop.");
  }

  observers_.reset();
  Rmat_.reset();
  qcflags_.clear();
  gradFG_.reset();
  yobs_.reset();

  obspaces_.appendObs(appendConfig);

  yobs_.reset(new Observations_(obspaces_, "ObsValue"));
  observers_.reset(new Observers_(obspaces_, conf_));
  Rmat_.reset(new ObsErrors_(eckit::LocalConfiguration(conf_, "observers"), obspaces_));
  for (size_t jj = 0; jj < obspaces_.size(); ++jj) {
    ObsDataInt_ qcflags(obspaces_[jj], obspaces_[jj].obsvariables());
    qcflags_.push_back(qcflags);
  }

  Log::trace() << "CostJo::appendObs done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJO_H_
