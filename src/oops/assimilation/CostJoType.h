/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_COSTJOTYPE_H_
#define OOPS_ASSIMILATION_COSTJOTYPE_H_

#include <memory>
#include <string>
#include <vector>

#include <boost/noncopyable.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/GetValuePost.h"
#include "oops/base/ObsErrorBase.h"
#include "oops/base/Observer.h"
#include "oops/base/PostBase.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/ObsAuxControl.h"
#include "oops/interface/ObsSpace.h"
#include "oops/interface/ObsVector.h"
#include "oops/interface/State.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Jo Cost Function
/*!
 * The CostJoType class encapsulates the Jo term of the cost function for a single obs type.
 */

template<typename MODEL, typename OBS> class CostJoType : private boost::noncopyable {
  typedef Geometry<MODEL>               Geometry_;
  typedef ObsAuxControl<OBS>            ObsAuxCtrl_;
  typedef ObsErrorBase<OBS>             ObsError_;
  typedef Observer<MODEL, OBS>          Observer_;
  typedef ObsSpace<OBS>                 ObsSpace_;
  typedef ObsVector<OBS>                ObsVector_;
  typedef PostBase<State<MODEL>>        PostBase_;
  typedef State<MODEL>                  State_;
  typedef std::shared_ptr<GetValuePost<MODEL, OBS>> GetValuePtr_;

 public:
  /// Construct \f$ J_o\f$ from \f$ R\f$ and \f$ y_{obs}\f$.
  CostJoType(ObsSpace_ &, const eckit::Configuration &);

  /// Destructor
  virtual ~CostJoType() {}

  /// Initialize \f$ J_o\f$ before starting the integration of the model.
  GetValuePtr_ initialize(const Geometry_ &, const ObsAuxCtrl_ &, const eckit::Configuration &);

  /// Finalize \f$ J_o\f$ after the integration of the model.
  void finalize(ObsVector_ &, ObsVector_ &, ObsVector_ &, ObsVector_ &);

  /// Multiply by \f$ R\f$ and \f$ R^{-1}\f$.
  void multiplyR(ObsVector_ &) const;
  void inverseMultiplyR(ObsVector_ &) const;

  /// Print Jo
  double printJo(const ObsVector_ &, const ObsVector_ &) const;

 private:
  eckit::LocalConfiguration obsconf_;
  ObsSpace_ & obspace_;
  ObsVector_ obserr_;  // Obs errors (should come from R matrix)

  std::unique_ptr<ObsError_> Rmat_;

  /// Configuration for current initialize/finalize pair
  std::unique_ptr<eckit::Configuration> currentConf_;
  int iter_;

  /// Used for computing H(x) and running QC filters
  Observer_ observer_;
};

// =============================================================================

template<typename MODEL, typename OBS>
CostJoType<MODEL, OBS>::CostJoType(ObsSpace_ & obspace, const eckit::Configuration & conf)
  : obsconf_(conf), obspace_(obspace),
    obserr_(obspace_, "ObsError"),
    Rmat_(), currentConf_(), iter_(0), observer_(obspace_, obsconf_)
{
  Log::trace() << "CostJoType::CostJoType done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
std::shared_ptr<GetValuePost<MODEL, OBS>>
CostJoType<MODEL, OBS>::initialize(const Geometry_ & geom, const ObsAuxCtrl_ & ybias,
                                   const eckit::Configuration & conf) {
  Log::trace() << "CostJoType::initialize start" << std::endl;

  currentConf_.reset(new eckit::LocalConfiguration(conf));
  iter_ = currentConf_->getInt("iteration");

  GetValuePtr_ getvals = observer_.initialize(geom, ybias, obserr_, iter_);

  Log::trace() << "CostJoType::initialize done" << std::endl;
  return getvals;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJoType<MODEL, OBS>::finalize(ObsVector_ & yobs, ObsVector_ & yeqv,
                                      ObsVector_ & ydep, ObsVector_ & grad) {
  Log::trace() << "CostJoType::finalize start" << std::endl;

// Get simulated obs values from observer
  observer_.finalize(yeqv);

// Sace current simulated obs, QC flags and obs error
  const std::string obsname = "hofx" + std::to_string(iter_);
  yeqv.save(obsname);

  const std::string errname = "EffectiveError" + std::to_string(iter_);
  obserr_.save(errname);
  obserr_.save("EffectiveError");  // Obs error covariance is looking for that for now

// Set observation error covariance
  const eckit::LocalConfiguration rconf(obsconf_, "obs error");
  Rmat_ = ObsErrorFactory<OBS>::create(rconf, obspace_);

// Perturb observations according to obs error statistics
  bool obspert = currentConf_->getBool("obs perturbations", false);
  if (obspert) {
    ObsVector_ ypert(obspace_);
    Rmat_->randomize(ypert);
    yobs += ypert;
    Log::info() << "Perturbed observations: " << yobs << std::endl;
  }

// Compute departures and Jo gradient
  ydep = yeqv;
  ydep -= yobs;
  grad = ydep;
  Rmat_->inverseMultiply(grad);

// Save departures for diagnostics if required
  if (currentConf_->has("diagnostics.departures")) {
    const std::string depname = currentConf_->getString("diagnostics.departures");
    ydep.save(depname);
  }

  currentConf_.reset();
  Log::trace() << "CostJoType::finalize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
double CostJoType<MODEL, OBS>::printJo(const ObsVector_ & dy, const ObsVector_ & grad) const {
  Log::trace() << "CostJoType::printJo start" << std::endl;

  double zjo = 0.0;
  const unsigned nobs = grad.nobs();
  Log::test() << "CostJo   : Nonlinear Jo(" << obspace_.obsname() << ") = ";

  if (nobs > 0) {
    zjo = 0.5 * dot_product(dy, grad);
    const double err = Rmat_->getRMSE();
    Log::test() << zjo << ", nobs = " << nobs << ", Jo/n = " << zjo/nobs << ", err = " << err;
  } else {
    Log::test() << zjo << " --- No Observations";
  }

  if (obsconf_.getBool("monitoring only", false)) {
    Log::test() << " (Monitoring only)";
    zjo = 0.0;
  }
  Log::test() << std::endl;

  Log::trace() << "CostJoType::printJo done" << std::endl;
  return zjo;
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJoType<MODEL, OBS>::multiplyR(ObsVector_ & dy) const {
  Rmat_->multiply(dy);
}

// -----------------------------------------------------------------------------

template<typename MODEL, typename OBS>
void CostJoType<MODEL, OBS>::inverseMultiplyR(ObsVector_ & dy) const {
  Rmat_->inverseMultiply(dy);
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_COSTJOTYPE_H_
