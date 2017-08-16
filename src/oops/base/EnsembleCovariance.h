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

#include "eckit/config/LocalConfiguration.h"
#include "util/Logger.h"
#include "oops/assimilation/GMRESR.h"
#include "oops/base/Ensemble.h"
#include "oops/base/EnsemblesCollection.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/Localization.h"
#include "oops/interface/State.h"
#include "oops/interface/Variables.h"
#include "util/DateTime.h"
#include "util/abor1_cpp.h"

namespace oops {

/// Generic ensemble based model space error covariance.

// -----------------------------------------------------------------------------
template <typename MODEL>
class EnsembleCovariance : public ModelSpaceCovarianceBase<MODEL> {
  typedef Ensemble<MODEL>                       Ensemble_;
  typedef boost::shared_ptr<Ensemble<MODEL> >   EnsemblePtr_;
  typedef EnsemblesCollection<MODEL>            EnsemblesCollection_;
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef Localization<MODEL>        Localization_;
  typedef State<MODEL>               State_;
  typedef Variables<MODEL>           Variables_;

 public:
  EnsembleCovariance(const Geometry_ &, const Variables_ &, const eckit::Configuration &, const State_ &);
  ~EnsembleCovariance();

  void linearize(const State_ &, const Geometry_ &) override;

  void multiply(const Increment_ &, Increment_ &) const override;
  void inverseMultiply(const Increment_ &, Increment_ &) const override;

  void randomize(Increment_ &) const override;

 private:
  const eckit::LocalConfiguration config_;
  const util::DateTime time_;
  boost::scoped_ptr<Localization_> loc_;
};

// =============================================================================

/// Constructor, destructor
// -----------------------------------------------------------------------------
template<typename MODEL>
EnsembleCovariance<MODEL>::EnsembleCovariance(const Geometry_ &, const Variables_ &,
                                              const eckit::Configuration & config, const State_ &)
  : config_(config), time_(config_.getString("date")), loc_()
{
  Log::trace() << "EnsembleCovariance created." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
EnsembleCovariance<MODEL>::~EnsembleCovariance() {
  Log::trace() << "EnsembleCovariance destructed." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void EnsembleCovariance<MODEL>::linearize(const State_ & xb,
                                          const Geometry_ & resol) {
  // Compute the ensemble of perturbations at time of xb.
  ASSERT(xb.validTime() == time_);
  EnsemblePtr_ ens_k(new Ensemble_(xb.validTime(), config_));
  ens_k->linearize(xb, resol);
  EnsemblesCollection_::getInstance().put(xb.validTime(), ens_k);

  const eckit::LocalConfiguration conf(config_, "localization");
  loc_.reset(new Localization_(xb, conf));
  Log::trace() << "EnsembleCovariance linearized." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void EnsembleCovariance<MODEL>::multiply(const Increment_ & dxi,
                                         Increment_ & dxo) const {
  EnsemblePtr_ e_k2 = EnsemblesCollection_::getInstance()[dxi.validTime()];
  EnsemblePtr_ e_k1 = EnsemblesCollection_::getInstance()[dxo.validTime()];

  // Apply to dxi the block B_k1k2 of the B matrix which contains the
  // covariance of perturbations at time k1 with covariance at time k2
  dxo.zero();
  for (unsigned int m = 0; m < e_k1->size(); ++m) {
    Increment_ dx(dxi);
    dx.schur_product_with((*e_k2)[m]);
    loc_->multiply(dx);
    dx.schur_product_with((*e_k1)[m]);

    // We can't use '+=' , dxo and dx aren't at the same time (QG)
    dxo.axpy(1.0, dx, false);
  }
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void EnsembleCovariance<MODEL>::inverseMultiply(const Increment_ & dxi,
                                                Increment_ & dxo) const {
  IdentityMatrix<Increment_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void EnsembleCovariance<MODEL>::randomize(Increment_ &) const {
  ABORT("EnsembleCovariance::randomize: Would it make sense?");
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_ENSEMBLECOVARIANCE_H_
