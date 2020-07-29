/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_HYBRIDCOVARIANCE_H_
#define OOPS_BASE_HYBRIDCOVARIANCE_H_

#include <memory>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/GMRESR.h"
#include "oops/base/EnsembleCovariance.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace oops {

/// Generic hybrid static-ensemble model space error covariance.

// -----------------------------------------------------------------------------
template <typename MODEL>
class HybridCovariance : public ModelSpaceCovarianceBase<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;

 public:
  HybridCovariance(const Geometry_ &, const Variables &,
                   const eckit::Configuration &, const State_ &, const State_ &);
  ~HybridCovariance();

 private:
  void doRandomize(Increment_ &) const override;
  void doMultiply(const Increment_ &, Increment_ &) const override;
  void doInverseMultiply(const Increment_ &, Increment_ &) const override;

  std::unique_ptr< ModelSpaceCovarianceBase<MODEL> > static_;
  std::unique_ptr< EnsembleCovariance<MODEL> >  ensemble_;
  double ensWeight_;
  double staWeight_;
};

// =============================================================================

/// Constructor, destructor
// -----------------------------------------------------------------------------
template<typename MODEL>
HybridCovariance<MODEL>::HybridCovariance(const Geometry_ & resol, const Variables & vars,
                                          const eckit::Configuration & config,
                                          const State_ & xb, const State_ & fg)
  : ModelSpaceCovarianceBase<MODEL>(xb, fg, resol, config),
    static_(CovarianceFactory<MODEL>::create(
              eckit::LocalConfiguration(config, "static"), resol, vars, xb, fg))
{
  const eckit::LocalConfiguration ensConf(config, "ensemble");
  ensemble_.reset(new EnsembleCovariance<MODEL>(resol,  vars, ensConf, xb, fg));

  ensWeight_ = config.getDouble("ensemble weight");
  ASSERT(ensWeight_ > 0.0);
  staWeight_ = config.getDouble("static weight");
  ASSERT(staWeight_ > 0.0);
  Log::trace() << "HybridCovariance created." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
HybridCovariance<MODEL>::~HybridCovariance() {
  Log::trace() << "HybridCovariance destructed" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void HybridCovariance<MODEL>::doMultiply(const Increment_ & dxi, Increment_ & dxo) const {
  static_->multiply(dxi, dxo);
  dxo *= staWeight_;
  Increment_ tmp(dxo);
  ensemble_->multiply(dxi, tmp);
  dxo.axpy(ensWeight_, tmp);
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void HybridCovariance<MODEL>::doInverseMultiply(const Increment_ & dxi, Increment_ & dxo) const {
  IdentityMatrix<Increment_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void HybridCovariance<MODEL>::doRandomize(Increment_ &) const {
  ABORT("HybridCovariance::doRandomize: Would it make sense?");
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_HYBRIDCOVARIANCE_H_
