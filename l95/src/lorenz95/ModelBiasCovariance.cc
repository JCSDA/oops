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
#include <random>
#include <string>

#include "util/Logger.h"
#include "lorenz95/ModelBiasCorrection.h"
#include "eckit/config/Configuration.h"

using oops::Log;

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
    Log::info() << "ModelBiasCovariance variance = " << variance_ << std::endl;
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
  static std::mt19937 generator(3);
  static std::normal_distribution<double> distribution(0.0, stdev);
  dx.bias() = distribution(generator);
}
// -----------------------------------------------------------------------------
void ModelBiasCovariance::print(std::ostream & os) const {
  os << "ModelBiasCovariance::print not implemented";
}
// -----------------------------------------------------------------------------
}  // namespace lorenz95
