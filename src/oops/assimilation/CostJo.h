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
#include "oops/assimilation/CostJoType.h"
#include "oops/assimilation/CostTermBase.h"
#include "oops/base/Departures.h"
#include "oops/base/GetValuePosts.h"
#include "oops/base/ObsErrors.h"
#include "oops/base/Observations.h"
#include "oops/base/ObserversTLAD.h"
#include "oops/base/ObsSpaces.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/PostProcessorTLAD.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
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
  typedef CostJoType<MODEL, OBS>        JoType_;
  typedef Departures<OBS>               Departures_;
  typedef Observations<OBS>             Observations_;
  typedef Geometry<MODEL>               Geometry_;
  typedef GetValuePosts<MODEL, OBS>     GetValuePosts_;
  typedef State<MODEL>                  State_;
  typedef Increment<MODEL>              Increment_;
  typedef ObsErrors<OBS>                ObsErrors_;
  typedef ObsSpaces<OBS>                ObsSpaces_;
  typedef ObserversTLAD<MODEL, OBS>     ObserversTLAD_;
  typedef PostProcessor<State_>         PostProc_;
  typedef PostProcessorTLAD<MODEL>      PostProcTLAD_;

 public:
  /// Construct \f$ J_o\f$ from \f$ R\f$ and \f$ y_{obs}\f$.
  CostJo(const eckit::Configuration &, const eckit::mpi::Comm &,
         const util::DateTime &, const util::DateTime &,
         const eckit::mpi::Comm & ctime = oops::mpi::myself());

  /// Destructor
  virtual ~CostJo() {}

  /// Initialize \f$ J_o\f$ before starting the integration of the model.
  void initialize(const CtrlVar_ &, const eckit::Configuration &, PostProc_ &) override;
  void initializeTraj(const CtrlVar_ &, const Geometry_ &,
                      const eckit::Configuration &, PostProcTLAD_ &) override;
  /// Finalize \f$ J_o\f$ after the integration of the model.
  double finalize() override;
  void finalizeTraj() override;

  /// Initialize \f$ J_o\f$ before starting the TL run.
  void setupTL(const CtrlInc_ &, PostProcTLAD_ &) const override;

  /// Initialize \f$ J_o\f$ before starting the AD run.
  void setupAD(std::shared_ptr<const GeneralizedDepartures>,
               CtrlInc_ &, PostProcTLAD_ &) const override;

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

  /// Accessor...
  const ObsSpaces_ & obspaces() const {return obspaces_;}

 private:
  const eckit::LocalConfiguration obsconf_;
  ObsSpaces_ obspaces_;
  Observations_ yobs_;

  std::vector<std::unique_ptr<JoType_>> jos_;

  /// Jo Gradient at first guess : \f$ R^{-1} (H(x_{fg})-y_{obs}) \f$.
  std::unique_ptr<Departures_> gradFG_;

  /// Linearized observation operators.
  std::shared_ptr<ObserversTLAD_> pobstlad_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostJo<MODEL, OBS>::CostJo(const eckit::Configuration & joConf, const eckit::mpi::Comm & comm,
                           const util::DateTime & winbgn, const util::DateTime & winend,
                           const eckit::mpi::Comm & ctime)
  : obsconf_(joConf), obspaces_(obsconf_, comm, winbgn, winend, ctime),
    yobs_(obspaces_, "ObsValue"), jos_(), gradFG_()
{
  Log::trace() << "CostJo::CostJo start" << std::endl;
  jos_.reserve(obspaces_.size());
  std::vector<eckit::LocalConfiguration> confs(obsconf_.getSubConfigurations());
  for (size_t jj = 0; jj < obspaces_.size(); ++jj) {
    jos_.emplace_back(new JoType_(obspaces_[jj], confs[jj]));
  }
  Log::trace() << "CostJo::CostJo done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJo<MODEL, OBS>::initialize(const CtrlVar_ & xx, const eckit::Configuration & conf,
                                    PostProc_ & pp) {
  Log::trace() << "CostJo::initialize start" << std::endl;
  gradFG_.reset();

  std::shared_ptr<GetValuePosts_> getvals(new GetValuePosts_());
  for (size_t jj = 0; jj < jos_.size(); ++jj) {
    getvals->append(jos_[jj]->initialize(xx.state().geometry(), xx.obsVar()[jj], conf));
  }
  pp.enrollProcessor(getvals);

  Log::trace() << "CostJo::initialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
double CostJo<MODEL, OBS>::finalize() {
  Log::trace() << "CostJo::finalize start" << std::endl;

  // Obs, simulated obs and departures (held here for nice prints and diagnostics)
  Observations_ yeqv(obspaces_);
  Departures_ ydep(obspaces_);

  // Gradient at first guess (to define inner loop rhs)
  gradFG_.reset(new Departures_(obspaces_));

  // Actual Jo computations
  for (size_t jj = 0; jj < jos_.size(); ++jj) {
    jos_[jj]->finalize(yobs_[jj], yeqv[jj], ydep[jj], (*gradFG_)[jj]);
  }

  // Print diagnostics
  Log::info() << "Jo Observations:" << std::endl << yobs_
              << "End Jo Observations" << std::endl;

  Log::info() << "Jo Observations Equivalent:" << std::endl << yeqv
              << "End Jo Observations Equivalent" << std::endl;

  Log::info() << "Jo Bias Corrected Departures:" << std::endl << ydep
          << "End Jo Bias Corrected Departures" << std::endl;

  // Print Jo table
  double zjo = 0.0;
  for (size_t jj = 0; jj < jos_.size(); ++jj) {
    zjo += jos_[jj]->printJo(ydep[jj], (*gradFG_)[jj]);
  }

  Log::trace() << "CostJo::finalize done" << std::endl;
  return zjo;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJo<MODEL, OBS>::initializeTraj(const CtrlVar_ & xx, const Geometry_ &,
                                        const eckit::Configuration & conf,
                                        PostProcTLAD_ & pptraj) {
  Log::trace() << "CostJo::initializeTraj start" << std::endl;
  pobstlad_.reset(new ObserversTLAD_(obsconf_, obspaces_, xx.obsVar()));
  pptraj.enrollProcessor(pobstlad_);
  Log::trace() << "CostJo::initializeTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJo<MODEL, OBS>::finalizeTraj() {
  Log::trace() << "CostJo::finalizeTraj done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJo<MODEL, OBS>::setupTL(const CtrlInc_ & dx, PostProcTLAD_ & pptl) const {
  Log::trace() << "CostJo::setupTL start" << std::endl;
  ASSERT(pobstlad_);
  pobstlad_->setupTL(dx.obsVar());
  pptl.enrollProcessor(pobstlad_);
  Log::trace() << "CostJo::setupTL done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJo<MODEL, OBS>::setupAD(std::shared_ptr<const GeneralizedDepartures> pv,
                                 CtrlInc_ & dx, PostProcTLAD_ & ppad) const {
  Log::trace() << "CostJo::setupAD start" << std::endl;
  ASSERT(pobstlad_);
  std::shared_ptr<const Departures_> dy = std::dynamic_pointer_cast<const Departures_>(pv);
  pobstlad_->setupAD(dy, dx.obsVar());
  ppad.enrollProcessor(pobstlad_);
  Log::trace() << "CostJo::setupAD done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
std::unique_ptr<GeneralizedDepartures>
CostJo<MODEL, OBS>::multiplyCovar(const GeneralizedDepartures & v1) const {
  Log::trace() << "CostJo::multiplyCovar start" << std::endl;
  std::unique_ptr<Departures_> y1(new Departures_(dynamic_cast<const Departures_ &>(v1)));
  for (size_t jj = 0; jj < jos_.size(); ++jj) {
    jos_[jj]->multiplyR((*y1)[jj]);
  }
  return std::move(y1);
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
std::unique_ptr<GeneralizedDepartures>
CostJo<MODEL, OBS>::multiplyCoInv(const GeneralizedDepartures & v1) const {
  Log::trace() << "CostJo::multiplyCoInv start" << std::endl;
  std::unique_ptr<Departures_> y1(new Departures_(dynamic_cast<const Departures_ &>(v1)));
  for (size_t jj = 0; jj < jos_.size(); ++jj) {
    jos_[jj]->inverseMultiplyR((*y1)[jj]);
  }
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
void CostJo<MODEL, OBS>::resetLinearization() {
  Log::trace() << "CostJo::resetLinearization start" << std::endl;
  pobstlad_.reset();
  Log::trace() << "CostJo::resetLinearization done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJO_H_
