/*
 * (C) Copyright 2025- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "oops/base/ModelSpaceCovarianceBase.h"
#include "oops/coupled/DataSetCoupled.h"
#include "oops/coupled/GeometryCoupled.h"
#include "oops/coupled/IncrementCoupled.h"
#include "oops/coupled/StateCoupled.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace oops {
template <typename MODEL1, typename MODEL2> struct TraitCoupled;

/// \brief Block diagonal covariance for coupled model, using CovarianceFactory.
/// \details This class implements a block diagonal covariance for coupled models,
/// where each block corresponds to the covariance for one of the component models.
/// The component model covariances are created using the CovarianceFactory for each model.
template<typename MODEL1, typename MODEL2>
class BlockDiagonalCovarianceCoupled :
    public ModelSpaceCovarianceBase<TraitCoupled<MODEL1, MODEL2>> {
  typedef TraitCoupled<MODEL1, MODEL2> COUPLED;
  typedef Geometry<COUPLED>                Geometry_;
  typedef Increment4D<COUPLED>             Increment4D_;
  typedef State4D<COUPLED>                 State4D_;
  typedef ModelSpaceCovarianceBase<MODEL1> ErrorCovariance1_;
  typedef ModelSpaceCovarianceBase<MODEL2> ErrorCovariance2_;
 public:
  static const std::string classname() {return "BlockDiagonalCovarianceCoupled";}

  BlockDiagonalCovarianceCoupled(const Geometry_ &, const oops::Variables &,
                         const eckit::Configuration &, const State4D_ &, const State4D_ &);
  ~BlockDiagonalCovarianceCoupled() = default;

 private:
  // Implementation of pure virtual methods required by ModelSpaceCovarianceBase
  void doMultiply(const Increment4D_ &, Increment4D_ &) const override;
  void doInverseMultiply(const Increment4D_ &, Increment4D_ &) const override;
  void doRandomize(Increment4D_ &) const override;

  // Component model covariances
  std::unique_ptr<ErrorCovariance1_> cov1_;
  std::unique_ptr<ErrorCovariance2_> cov2_;
};

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
BlockDiagonalCovarianceCoupled<MODEL1, MODEL2>::BlockDiagonalCovarianceCoupled(
    const Geometry_ & geom,
    const oops::Variables & vars,
    const eckit::Configuration & config,
    const State4D_ & xb,
    const State4D_ & fg):
    ModelSpaceCovarianceBase<TraitCoupled<MODEL1, MODEL2>>(geom, config, xb, fg) {
  Log::trace() << "BlockDiagonalCovarianceCoupled::ctor starting" << std::endl;
  const eckit::LocalConfiguration conf1(config, MODEL1::name());
  const eckit::LocalConfiguration conf2(config, MODEL2::name());
  const std::vector<Variables> splitvars = splitVariables(vars, geom.geometry().variables());

  // Create component State4D objects using specialized constructors (zero-copy)
  State4D<MODEL1> xb1 = share_state1(xb);
  State4D<MODEL2> xb2 = share_state2(xb);
  State4D<MODEL1> fg1 = share_state1(fg);
  State4D<MODEL2> fg2 = share_state2(fg);
  cov1_ = std::unique_ptr<ErrorCovariance1_>(CovarianceFactory<MODEL1>::create(
               geom.geometry().geometry1(), splitvars[0], conf1, xb1, fg1));
  cov2_ = std::unique_ptr<ErrorCovariance2_>(CovarianceFactory<MODEL2>::create(
               geom.geometry().geometry2(), splitvars[1], conf2, xb2, fg2));
  oops::Log::trace() << "BlockDiagonalCovarianceCoupled::ctor done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void BlockDiagonalCovarianceCoupled<MODEL1, MODEL2>::doMultiply(
    const Increment4D_ & dx, Increment4D_ & dy) const {
  Log::trace() << "BlockDiagonalCovarianceCoupled::doMultiply starting" << std::endl;
  // Create component increments using specialized constructors (no copies made)
  Increment4D<MODEL1> dx1 = share_increment1(dx);
  Increment4D<MODEL2> dx2 = share_increment2(dx);
  Increment4D<MODEL1> dy1 = share_increment1(dy);
  Increment4D<MODEL2> dy2 = share_increment2(dy);
  cov1_->multiply(dx1, dy1);
  cov2_->multiply(dx2, dy2);
  Log::trace() << "BlockDiagonalCovarianceCoupled::doMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void BlockDiagonalCovarianceCoupled<MODEL1, MODEL2>::doInverseMultiply(
    const Increment4D_ & dx, Increment4D_ & dy) const {
  Log::trace() << "BlockDiagonalCovarianceCoupled::doInverseMultiply starting" << std::endl;
  // Create component increments using specialized constructors (no copies made)
  Increment4D<MODEL1> dx1 = share_increment1(dx);
  Increment4D<MODEL2> dx2 = share_increment2(dx);
  Increment4D<MODEL1> dy1 = share_increment1(dy);
  Increment4D<MODEL2> dy2 = share_increment2(dy);
  cov1_->inverseMultiply(dx1, dy1);
  cov2_->inverseMultiply(dx2, dy2);
  Log::trace() << "BlockDiagonalCovarianceCoupled::doInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void BlockDiagonalCovarianceCoupled<MODEL1, MODEL2>::doRandomize(Increment4D_ & dx) const {
  Log::trace() << "BlockDiagonalCovarianceCoupled::doRandomize starting" << std::endl;
  // Create component increments using specialized constructors (no copies made)
  Increment4D<MODEL1> dx1 = share_increment1(dx);
  Increment4D<MODEL2> dx2 = share_increment2(dx);
  cov1_->randomize(dx1);
  cov2_->randomize(dx2);
  Log::trace() << "BlockDiagonalCovarianceCoupled::doRandomize done" << std::endl;
}

}  // namespace oops
