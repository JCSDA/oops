/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_ENSEMBLECOVARIANCE_H_
#define OOPS_BASE_ENSEMBLECOVARIANCE_H_

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/system/ResourceUsage.h"

#include "oops/assimilation/GMRESR.h"
#include "oops/base/Geometry.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/Increment.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/Localization.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/util/Logger.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Timer.h"

namespace oops {

/// Parameters describing generic ensemble covariances.
template <typename MODEL>
class EnsembleCovarianceParameters : public ModelSpaceCovarianceParametersBase<MODEL> {
  OOPS_CONCRETE_PARAMETERS(EnsembleCovarianceParameters,
                           ModelSpaceCovarianceParametersBase<MODEL>)
 public:
  /// Parameters for ensemble of increments used in the covariances.
  IncrementEnsembleFromStatesParameters<MODEL> ensemble{this};
  oops::OptionalParameter<eckit::LocalConfiguration> localization{"localization",
                         "localization applied to ensemble covariances", this};
};

/// Generic ensemble based model space error covariance.

// -----------------------------------------------------------------------------
template <typename MODEL>
class EnsembleCovariance : public ModelSpaceCovarianceBase<MODEL>,
                           private util::ObjectCounter<EnsembleCovariance<MODEL>> {
  typedef Geometry<MODEL>                           Geometry_;
  typedef Increment<MODEL>                          Increment_;
  typedef Localization<MODEL>                       Localization_;
  typedef State<MODEL>                              State_;
  typedef IncrementEnsemble<MODEL>                  Ensemble_;
  typedef std::shared_ptr<IncrementEnsemble<MODEL>> EnsemblePtr_;

 public:
  typedef EnsembleCovarianceParameters<MODEL> Parameters_;

  static const std::string classname() {return "oops::EnsembleCovariance";}

  EnsembleCovariance(const Geometry_ &, const Variables &,
                     const Parameters_ &, const State_ &, const State_ &);
  ~EnsembleCovariance();

 private:
  void doRandomize(Increment_ &) const override;
  void doMultiply(const Increment_ &, Increment_ &) const override;
  void doInverseMultiply(const Increment_ &, Increment_ &) const override;

  EnsemblePtr_ ens_;
  std::unique_ptr<Localization_> loc_;
  int seed_ = 7;  // For reproducibility
};

// =============================================================================

/// Constructor, destructor
// -----------------------------------------------------------------------------
template<typename MODEL>
EnsembleCovariance<MODEL>::EnsembleCovariance(const Geometry_ & resol, const Variables & vars,
                                              const Parameters_ & params,
                                              const State_ & xb, const State_ & fg)
  : ModelSpaceCovarianceBase<MODEL>(resol, params, xb, fg), ens_(), loc_()
{
  Log::trace() << "EnsembleCovariance::EnsembleCovariance start" << std::endl;
  util::Timer timer("oops::Covariance", "EnsembleCovariance");
  size_t init = eckit::system::ResourceUsage().maxResidentSetSize();
  ens_.reset(new Ensemble_(params.ensemble, xb, fg, resol, vars));
  if (ens_->size() < 2) {
    throw eckit::BadParameter("Not enough ensemble members provided for ensemble "
                              "covariances (at least 2 required)", Here());
  }
  if (params.localization.value() != boost::none) {
    eckit::LocalConfiguration conf = *params.localization.value();
    conf.set("date", xb.validTime().toString());
    loc_.reset(new Localization_(resol, xb.variables(), conf));
  }
  size_t current = eckit::system::ResourceUsage().maxResidentSetSize();
  this->setObjectSize(current - init);
  Log::trace() << "EnsembleCovariance::EnsembleCovariance done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
EnsembleCovariance<MODEL>::~EnsembleCovariance() {
  Log::trace() << "EnsembleCovariance destructed." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void EnsembleCovariance<MODEL>::doRandomize(Increment_ & dx) const {
  dx.zero();
  if (loc_) {
    // Localized covariance matrix
    for (unsigned int ie = 0; ie < ens_->size(); ++ie) {
      Increment_ tmp(dx);
      loc_->randomize(tmp);
      tmp.schur_product_with((*ens_)[ie]);
      dx.axpy(1.0, tmp, false);
    }
  } else {
    // Raw covariance matrix
    util::NormalDistribution<double> normalDist(ens_->size(), 0.0, 1.0, seed_);
    for (unsigned int ie = 0; ie < ens_->size(); ++ie) {
      dx.axpy(normalDist[ie], (*ens_)[ie]);
    }
  }
  dx *= 1.0/sqrt(static_cast<double>(ens_->size()-1));
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void EnsembleCovariance<MODEL>::doMultiply(const Increment_ & dxi, Increment_ & dxo) const {
  dxo.zero();
  for (unsigned int ie = 0; ie < ens_->size(); ++ie) {
    if (loc_) {
      // Localized covariance matrix
      Increment_ dx(dxi);
      dx.schur_product_with((*ens_)[ie]);
      loc_->multiply(dx);
      dx.schur_product_with((*ens_)[ie]);
      dxo.axpy(1.0, dx, false);
    } else {
      // Raw covariance matrix
      double wgt = dxi.dot_product_with((*ens_)[ie]);
      dxo.axpy(wgt, (*ens_)[ie], false);
    }
  }
  const double rk = 1.0/static_cast<double>(ens_->size()-1);
  dxo *= rk;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void EnsembleCovariance<MODEL>::doInverseMultiply(const Increment_ & dxi, Increment_ & dxo) const {
  IdentityMatrix<Increment_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_ENSEMBLECOVARIANCE_H_
