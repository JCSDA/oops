/*
 * (C) Copyright 2019-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */


#ifndef OOPS_BASE_HYBRIDCOVARIANCE4D_H_
#define OOPS_BASE_HYBRIDCOVARIANCE4D_H_

#include <memory>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/GMRESR.h"
#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/EnsembleCovariance4D.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/ModelSpaceCovariance4DBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace oops {

/// Generic hybrid static-ensemble model space error covariance.

// -----------------------------------------------------------------------------
template <typename MODEL>
class HybridCovariance4D : public ModelSpaceCovariance4DBase<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef Increment4D<MODEL>         Increment4D_;
  typedef State4D<MODEL>             State4D_;

 public:
  HybridCovariance4D(const Geometry_ &, const Variables &,
                   const eckit::Configuration &, const State4D_ &, const State4D_ &);
  ~HybridCovariance4D();

 private:
  void doRandomize(Increment4D_ &) const override;
  void doMultiply(const Increment4D_ &, Increment4D_ &) const override;
  void doInverseMultiply(const Increment4D_ &, Increment4D_ &) const override;

  std::unique_ptr< ModelSpaceCovarianceBase<MODEL> > static_;
  std::unique_ptr< EnsembleCovariance4D<MODEL> >  ensemble_;
  double ensWeight_;
  double staWeight_;
};

// =============================================================================

/// Constructor, destructor
// -----------------------------------------------------------------------------
template<typename MODEL>
HybridCovariance4D<MODEL>::HybridCovariance4D(const Geometry_ & resol, const Variables & vars,
                                              const eckit::Configuration & config,
                                              const State4D_ & xb, const State4D_ & fg)
  : ModelSpaceCovariance4DBase<MODEL>(xb, fg, resol, config),
    static_(CovarianceFactory<MODEL>::create(
              eckit::LocalConfiguration(config, "static"), resol, vars, xb[0], fg[0]))
{
  const eckit::LocalConfiguration ensConf(config, "ensemble");
  ensemble_.reset(new EnsembleCovariance4D<MODEL>(resol,  vars, ensConf, xb, fg));

  ensWeight_ = config.getDouble("ensemble weight");
  ASSERT(ensWeight_ > 0.0);
  staWeight_ = config.getDouble("static weight");
  ASSERT(staWeight_ > 0.0);
  Log::trace() << "HybridCovariance4D created." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
HybridCovariance4D<MODEL>::~HybridCovariance4D() {
  Log::trace() << "HybridCovariance4D destructed" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void HybridCovariance4D<MODEL>::doMultiply(const Increment4D_ & dxi, Increment4D_ & dxo) const {
  for (int isub = dxi.first(); isub <= dxi.last(); ++isub) {
    static_->multiply(dxi[isub], dxo[isub]);
  }
  dxo *= staWeight_;
  Increment4D_ tmp(dxo);
  ensemble_->multiply(dxi, tmp);
  dxo.axpy(ensWeight_, tmp);
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void HybridCovariance4D<MODEL>::doInverseMultiply(const Increment4D_ & dxi, Increment4D_ & dxo)
const {
  IdentityMatrix<Increment4D_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void HybridCovariance4D<MODEL>::doRandomize(Increment4D_ &) const {
  ABORT("HybridCovariance4D::doRandomize: Would it make sense?");
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_HYBRIDCOVARIANCE4D_H_
