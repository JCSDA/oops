/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_ENSEMBLECOVARIANCE4D_H_
#define OOPS_BASE_ENSEMBLECOVARIANCE4D_H_

#include <memory>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/GMRESR.h"
#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/LocalizationBase.h"
#include "oops/base/ModelSpaceCovariance4DBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace oops {

/// Generic ensemble based model space error 4D covariance.

// -----------------------------------------------------------------------------
template <typename MODEL>
class EnsembleCovariance4D : public ModelSpaceCovariance4DBase<MODEL> {
  typedef Geometry<MODEL>                           Geometry_;
  typedef Increment<MODEL>                          Increment_;
  typedef Increment4D<MODEL>                        Increment4D_;
  typedef LocalizationBase<MODEL>                   Localization_;
  typedef State4D<MODEL>                            State4D_;
  typedef IncrementEnsemble<MODEL>                  Ensemble_;
  typedef std::shared_ptr<IncrementEnsemble<MODEL>> EnsemblePtr_;

 public:
  EnsembleCovariance4D(const Geometry_ &, const Variables &,
                       const eckit::Configuration &, const State4D_ &, const State4D_ &);
  ~EnsembleCovariance4D();

 private:
  void doRandomize(Increment4D_ &) const override;
  void doMultiply(const Increment4D_ &, Increment4D_ &) const override;
  void doInverseMultiply(const Increment4D_ &, Increment4D_ &) const override;

  EnsemblePtr_ ens_;
  std::unique_ptr<Localization_> loc_;
};

// =============================================================================

/// Constructor, destructor
// -----------------------------------------------------------------------------
template<typename MODEL>
EnsembleCovariance4D<MODEL>::EnsembleCovariance4D(const Geometry_ & resol, const Variables & vars,
                                                  const eckit::Configuration & conf,
                                                  const State4D_ & xb, const State4D_ & fg)
  : ModelSpaceCovariance4DBase<MODEL>(xb, fg, resol, conf), ens_(), loc_()
{
  Log::trace() << "EnsembleCovariance4D::EnsembleCovariance4D start" << std::endl;
  ens_.reset(new Ensemble_(conf, xb, fg, resol, vars));
  if (conf.has("localization")) {
    const eckit::LocalConfiguration confloc(conf, "localization");
    loc_ = LocalizationFactory<MODEL>::create(resol, ens_, confloc);
  }
  Log::trace() << "EnsembleCovariance4D::EnsembleCovariance4D done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
EnsembleCovariance4D<MODEL>::~EnsembleCovariance4D() {
  Log::trace() << "EnsembleCovariance4D destructed." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void EnsembleCovariance4D<MODEL>::doRandomize(Increment4D_ &) const {
  ABORT("EnsembleCovariance4D::doRandomize: Would it make sense?");
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void EnsembleCovariance4D<MODEL>::doMultiply(const Increment4D_ & dxi,
                                             Increment4D_ & dxo) const {
  dxo.zero();
  for (unsigned int ie = 0; ie < ens_->size(); ++ie) {
    if (loc_) {
      // Localized covariance matrix
      Increment4D_ dx(dxi);
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
  const double rk = 1.0/(static_cast<double>(ens_->size()) - 1.0);
  dxo *= rk;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void EnsembleCovariance4D<MODEL>::doInverseMultiply(const Increment4D_ & dxi,
                                                    Increment4D_ & dxo) const {
  IdentityMatrix<Increment4D_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_ENSEMBLECOVARIANCE4D_H_

