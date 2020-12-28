/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/ModelBiasCovariance.h"

#include <cmath>
#include <iostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "lorenz95/ModelBiasCorrection.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------
ModelBiasCovariance::ModelBiasCovariance(const eckit::Configuration & conf, const Resolution &)
  : conf_(conf), variance_(0.0), active_(false)
{
  if (conf_.has("standard_deviation")) {
    active_ = true;
    const double zz = conf_.getDouble("standard_deviation");
    variance_ = zz * zz;
    ASSERT(variance_ > 0.0);
    oops::Log::info() << "ModelBiasCovariance variance = " << variance_ << std::endl;
  }
}
// -----------------------------------------------------------------------------
void ModelBiasCovariance::multiply(const ModelBiasCorrection & dxin,
                                   ModelBiasCorrection & dxout) const {
  if (active_) {
    dxout = dxin;
    dxout *= variance_;
  } else {
    dxout.zero();
  }
}
// -----------------------------------------------------------------------------
void ModelBiasCovariance::inverseMultiply(const ModelBiasCorrection & dxin,
                                          ModelBiasCorrection & dxout) const {
  if (active_) {
    dxout = dxin;
    dxout *= 1.0 / variance_;
  } else {
    dxout.zero();
  }
}
// -----------------------------------------------------------------------------
void ModelBiasCovariance::randomize(ModelBiasCorrection & dx) const {
  const double stdev = std::sqrt(variance_);
  util::NormalDistribution<double> x(1, 0.0, stdev, 3);
  dx.bias() = x[0];
}
// -----------------------------------------------------------------------------
void ModelBiasCovariance::print(std::ostream & os) const {
  if (active_) {
    os << "ModelBiasCovariance: variance = " << variance_;
  } else {
    os << "ModelBiasCovariance not active" << std::endl;
  }
}
// -----------------------------------------------------------------------------
}  // namespace lorenz95
