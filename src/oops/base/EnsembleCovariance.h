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
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/GMRESR.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/IdentityMatrix.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/Localization.h"
#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace oops {

/// Generic ensemble based model space error covariance.

// -----------------------------------------------------------------------------
template <typename MODEL>
class EnsembleCovariance : public ModelSpaceCovarianceBase<MODEL> {
  typedef Geometry<MODEL>                           Geometry_;
  typedef Increment<MODEL>                          Increment_;
  typedef Localization<MODEL>                       Localization_;
  typedef State<MODEL>                              State_;
  typedef State4D<MODEL>                            State4D_;
  typedef IncrementEnsemble<MODEL>                  Ensemble_;
  typedef std::shared_ptr<IncrementEnsemble<MODEL>> EnsemblePtr_;

 public:
  EnsembleCovariance(const Geometry_ &, const Variables &,
                     const eckit::Configuration &, const State_ &, const State_ &);
  ~EnsembleCovariance();

 private:
  void doRandomize(Increment_ &) const override;
  void doMultiply(const Increment_ &, Increment_ &) const override;
  void doInverseMultiply(const Increment_ &, Increment_ &) const override;

  EnsemblePtr_ ens_;
  std::unique_ptr<Localization_> loc_;
};

// =============================================================================

/// Constructor, destructor
// -----------------------------------------------------------------------------
template<typename MODEL>
EnsembleCovariance<MODEL>::EnsembleCovariance(const Geometry_ & resol, const Variables &,
                                              const eckit::Configuration & conf,
                                              const State_ & xb, const State_ & fg)
  : ModelSpaceCovarianceBase<MODEL>(xb, fg, resol, conf), ens_(), loc_()
{
  Log::trace() << "EnsembleCovariance::EnsembleCovariance start" << std::endl;
  State4D_ xb4D(xb);
  State4D_ fg4D(fg);
  ens_.reset(new Ensemble_(conf, xb4D, fg4D, resol));
  if (conf.has("localization")) {
    const eckit::LocalConfiguration confloc(conf, "localization");
    loc_.reset(new Localization_(resol, ens_, confloc));
  }
  Log::trace() << "EnsembleCovariance::EnsembleCovariance done" << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
EnsembleCovariance<MODEL>::~EnsembleCovariance() {
  Log::trace() << "EnsembleCovariance destructed." << std::endl;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void EnsembleCovariance<MODEL>::doMultiply(const Increment_ & dxi,
                                           Increment_ & dxo) const {
  dxo.zero();
  for (unsigned int ie = 0; ie < ens_->size(); ++ie) {
    if (loc_) {
      // Localized covariance matrix
      Increment_ dx(dxi);
      dx.schur_product_with((*ens_)[ie][0]);
      loc_->multiply(dx);
      dx.schur_product_with((*ens_)[ie][0]);
      dxo.axpy(1.0, dx, false);
    } else {
      // Raw covariance matrix
      double wgt = dxi.dot_product_with((*ens_)[ie][0]);
      dxo.axpy(wgt, (*ens_)[ie][0], false);
    }
  }
  const double rk = 1.0/(static_cast<double>(ens_->size()) - 1.0);
  dxo *= rk;
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void EnsembleCovariance<MODEL>::doInverseMultiply(const Increment_ & dxi,
                                                  Increment_ & dxo) const {
  IdentityMatrix<Increment_> Id;
  dxo.zero();
  GMRESR(dxo, dxi, *this, Id, 10, 1.0e-3);
}
// -----------------------------------------------------------------------------
template<typename MODEL>
void EnsembleCovariance<MODEL>::doRandomize(Increment_ &) const {
  ABORT("EnsembleCovariance::doRandomize: Would it make sense?");
}
// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_ENSEMBLECOVARIANCE_H_
