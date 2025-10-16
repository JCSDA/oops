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

#include "oops/coupled/GeometryCoupled.h"
#include "oops/coupled/IncrementCoupled.h"
#include "oops/coupled/StateCoupled.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"


// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace oops {

/// \brief Block diagonal covariance for coupled model, using MODEL implementations.
/// \details This class implements a block diagonal covariance for coupled models,
/// where each block corresponds to the covariance for one of the component models.
/// The component model covariances are created using the Covariance class defined
/// in each model trait.
template<typename MODEL1, typename MODEL2>
class ErrorCovarianceCoupled : public util::Printable {
  typedef GeometryCoupled<MODEL1, MODEL2>  GeometryCoupled_;
  typedef IncrementCoupled<MODEL1, MODEL2> IncrementCoupled_;
  typedef StateCoupled<MODEL1, MODEL2>     StateCoupled_;
  typedef typename MODEL1::Covariance      ErrorCovariance1_;
  typedef typename MODEL2::Covariance      ErrorCovariance2_;
 public:
  static const std::string classname() {return "ErrorCovarianceCoupled";}

  ErrorCovarianceCoupled(const GeometryCoupled_ &, const oops::Variables &,
                         const eckit::Configuration &,
                         const StateCoupled_ &, const StateCoupled_ &);
  ~ErrorCovarianceCoupled() = default;

  void multiply(const IncrementCoupled_ &, IncrementCoupled_ &) const;
  void inverseMultiply(const IncrementCoupled_ &, IncrementCoupled_ &) const;
  void randomize(IncrementCoupled_ &) const;
 private:
  std::unique_ptr<ErrorCovariance1_> cov1_;
  std::unique_ptr<ErrorCovariance2_> cov2_;
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
ErrorCovarianceCoupled<MODEL1, MODEL2>::ErrorCovarianceCoupled(
    const GeometryCoupled_ & geom,
    const oops::Variables & vars,
    const eckit::Configuration & config,
    const StateCoupled_ & xb,
    const StateCoupled_ & fg) {
  const eckit::LocalConfiguration conf1(config, MODEL1::name());
  const eckit::LocalConfiguration conf2(config, MODEL2::name());
  cov1_.reset(new ErrorCovariance1_(geom.geometry1().geometry(),
                  vars, conf1, xb.state1().state(), fg.state1().state()));
  cov2_.reset(new ErrorCovariance2_(geom.geometry2().geometry(),
                  vars, conf2, xb.state2().state(), fg.state2().state()));
  Log::trace() << "ErrorCovarianceCoupled created" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void ErrorCovarianceCoupled<MODEL1, MODEL2>::multiply(
    const IncrementCoupled_ & dx, IncrementCoupled_ & dy) const {
  Log::trace() << "ErrorCovarianceCoupled::multiply starting" << std::endl;
  cov1_->multiply(dx.increment1().increment(), dy.increment1().increment());
  cov2_->multiply(dx.increment2().increment(), dy.increment2().increment());
  Log::trace() << "ErrorCovarianceCoupled::multiply done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void ErrorCovarianceCoupled<MODEL1, MODEL2>::inverseMultiply(
    const IncrementCoupled_ & dx, IncrementCoupled_ & dy) const {
  Log::trace() << "ErrorCovarianceCoupled::inverseMultiply starting" << std::endl;
  cov1_->inverseMultiply(dx.increment1().increment(), dy.increment1().increment());
  cov2_->inverseMultiply(dx.increment2().increment(), dy.increment2().increment());
  Log::trace() << "ErrorCovarianceCoupled::inverseMultiply done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void ErrorCovarianceCoupled<MODEL1, MODEL2>::randomize(IncrementCoupled_ & dx) const {
  Log::trace() << "ErrorCovarianceCoupled::randomize starting" << std::endl;
  cov1_->randomize(dx.increment1().increment());
  cov2_->randomize(dx.increment2().increment());
  Log::trace() << "ErrorCovarianceCoupled::randomize done" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename MODEL1, typename MODEL2>
void ErrorCovarianceCoupled<MODEL1, MODEL2>::print(std::ostream & os) const {
  os << "ErrorCovarianceCoupled" << std::endl;
}

}  // namespace oops
