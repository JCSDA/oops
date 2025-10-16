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
  typedef Increment4D<MODEL1>              Increment4D_1_;
  typedef Increment4D<MODEL2>              Increment4D_2_;
  typedef State4D<COUPLED>                 State4D_;
  typedef State4D<MODEL1>                  State4D_1_;
  typedef State4D<MODEL2>                  State4D_2_;
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

  // Create proper State4D objects for each component model
  State4D_1_ xb1(xb.geometry().geometry().geometry1(),
                 xb[0].state().state1().variables(), xb.validTimes(), xb.commTime());
  State4D_2_ xb2(xb.geometry().geometry().geometry2(),
                 xb[0].state().state2().variables(), xb.validTimes(), xb.commTime());

  // Create proper FG State4D objects for each component model
  State4D_1_ fg1(fg.geometry().geometry().geometry1(),
                 fg[0].state().state1().variables(), fg.validTimes(), fg.commTime());
  State4D_2_ fg2(fg.geometry().geometry().geometry2(),
                 fg[0].state().state2().variables(), fg.validTimes(), fg.commTime());

  // Copy data from coupled State4D objects to component State4D objects
  for (size_t jt = 0; jt < xb.size(); ++jt) {
    // Get the state for this time
    const auto& xbState = xb[jt];
    const auto& fgState = fg[jt];
    // Copy component model states
    xb1[jt] = xbState.state().state1();
    xb2[jt] = xbState.state().state2();
    fg1[jt] = fgState.state().state1();
    fg2[jt] = fgState.state().state2();
  }
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
  Increment4D_1_ dx1(dx[0].geometry().geometry().geometry1(),
            dx[0].increment().increment1().variables(), dx.validTimes(), dx.commTime());
  Increment4D_1_ dy1(dy[0].geometry().geometry().geometry1(),
            dy[0].increment().increment1().variables(), dx.validTimes(), dx.commTime());
  // Copy data from coupled increment to component increment for all timesteps
  for (size_t jt = 0; jt < dx.size(); ++jt) {
    dx1[jt] = dx[jt].increment().increment1();
    dy1[jt] = dy[jt].increment().increment1();
  }
  // Apply the covariance to all timesteps at once
  cov1_->multiply(dx1, dy1);
  // Copy results back to coupled increments
  for (size_t jt = 0; jt < dx.size(); ++jt) {
    dy[jt].increment().increment1() = dy1[jt];
  }
  Increment4D_2_ dx2(dx[0].geometry().geometry().geometry2(),
         dx[0].increment().increment2().variables(), dx.validTimes(), dx.commTime());
  Increment4D_2_ dy2(dy[0].geometry().geometry().geometry2(),
         dy[0].increment().increment2().variables(), dx.validTimes(), dx.commTime());
  // Copy data from coupled increment to component increment for all timesteps
  for (size_t jt = 0; jt < dx.size(); ++jt) {
    dx2[jt] = dx[jt].increment().increment2();
    dy2[jt] = dy[jt].increment().increment2();
  }

  // Apply the covariance to all timesteps at once
  cov2_->multiply(dx2, dy2);
  // Copy results back to coupled increments
  for (size_t jt = 0; jt < dx.size(); ++jt) {
    dy[jt].increment().increment2() = dy2[jt];
  }
  Log::trace() << "BlockDiagonalCovarianceCoupled::doMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void BlockDiagonalCovarianceCoupled<MODEL1, MODEL2>::doInverseMultiply(
    const Increment4D_ & dx, Increment4D_ & dy) const {
  Log::trace() << "BlockDiagonalCovarianceCoupled::doInverseMultiply starting" << std::endl;

  Increment4D_1_ dx1(dx[0].geometry().geometry().geometry1(),
          dx[0].increment().increment1().variables(), dx.validTimes(), dx.commTime());
  Increment4D_1_ dy1(dy[0].geometry().geometry().geometry1(),
          dy[0].increment().increment1().variables(), dx.validTimes(), dx.commTime());
  // Copy data from coupled increment to component increment for all timesteps
  for (size_t jt = 0; jt < dx.size(); ++jt) {
    dx1[jt] = dx[jt].increment().increment1();
  }
  // Apply the inverse covariance to all timesteps at once
  cov1_->inverseMultiply(dx1, dy1);
  // Copy results back to coupled increments
  for (size_t jt = 0; jt < dx.size(); ++jt) {
    dy[jt].increment().increment1() = dy1[jt];
  }
  Increment4D_2_ dx2(dx[0].geometry().geometry().geometry2(),
            dx[0].increment().increment2().variables(), dx.validTimes(), dx.commTime());
  Increment4D_2_ dy2(dy[0].geometry().geometry().geometry2(),
            dy[0].increment().increment2().variables(), dx.validTimes(), dx.commTime());
  // Copy data from coupled increment to component increment for all timesteps
  for (size_t jt = 0; jt < dx.size(); ++jt) {
    dx2[jt] = dx[jt].increment().increment2();
  }
  // Apply the inverse covariance to all timesteps at once
  cov2_->inverseMultiply(dx2, dy2);
  // Copy results back to coupled increments
  for (size_t jt = 0; jt < dx.size(); ++jt) {
    dy[jt].increment().increment2() = dy2[jt];
  }
  Log::trace() << "BlockDiagonalCovarianceCoupled::doInverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void BlockDiagonalCovarianceCoupled<MODEL1, MODEL2>::doRandomize(Increment4D_ & dx) const {
  Log::trace() << "BlockDiagonalCovarianceCoupled::doRandomize starting" << std::endl;
  Increment4D_1_ dx1(dx[0].geometry().geometry().geometry1(),
        dx[0].increment().increment1().variables(), dx.validTimes(), dx.commTime());
  // Randomize the entire component increment at once
  cov1_->randomize(dx1);
  // Copy results back to coupled increments
  for (size_t jt = 0; jt < dx.size(); ++jt) {
    dx[jt].increment().increment1() = dx1[jt];
  }
  Increment4D_2_ dx2(dx[0].geometry().geometry().geometry2(),
        dx[0].increment().increment2().variables(), dx.validTimes(), dx.commTime());
  // Randomize the entire component increment at once
  cov2_->randomize(dx2);
  // Copy results back to coupled increments
  for (size_t jt = 0; jt < dx.size(); ++jt) {
    dx[jt].increment().increment2() = dx2[jt];
  }
  Log::trace() << "BlockDiagonalCovarianceCoupled::doRandomize done" << std::endl;
}

}  // namespace oops
