/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/ModelBiasCorrection.h"

#include <iostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelBiasCovariance.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------
ModelBiasCorrection::ModelBiasCorrection(const Resolution &, const eckit::Configuration & conf)
  : bias_(0.0), active_(conf.has("standard_deviation"))
{
  if (active_) {
    oops::Log::trace() << "ModelBiasCorrection::ModelBiasCorrection created." << std::endl;
  }
}
// -----------------------------------------------------------------------------
ModelBiasCorrection::ModelBiasCorrection(const ModelBiasCorrection & other,
                                         const bool copy)
  : bias_(0.0), active_(other.active_)
{
  if (active_ && copy) bias_ = other.bias_;
}
// -----------------------------------------------------------------------------
ModelBiasCorrection::ModelBiasCorrection(const ModelBiasCorrection & other,
                                         const eckit::Configuration &)
  : bias_(0.0), active_(other.active_)
{
  if (active_) bias_ = other.bias_;
}
// -----------------------------------------------------------------------------
void ModelBiasCorrection::diff(const ModelBias & b1, const ModelBias & b2) {
  if (active_) bias_ = b1.bias() - b2.bias();
}
// -----------------------------------------------------------------------------
void ModelBiasCorrection::zero() {
  bias_ = 0.0;
}
// -----------------------------------------------------------------------------
ModelBiasCorrection & ModelBiasCorrection::operator=(const ModelBiasCorrection & rhs) {
  if (active_) {
    bias_ = rhs.bias_;
  } else {
    bias_ = 0.0;
  }
  return *this;
}
// -----------------------------------------------------------------------------
ModelBiasCorrection & ModelBiasCorrection::operator+=(const ModelBiasCorrection & rhs) {
  if (active_) bias_ += rhs.bias_;
  return *this;
}
// -----------------------------------------------------------------------------
ModelBiasCorrection & ModelBiasCorrection::operator-=(const ModelBiasCorrection & rhs) {
  if (active_) bias_ -= rhs.bias_;
  return *this;
}
// -----------------------------------------------------------------------------
ModelBiasCorrection & ModelBiasCorrection::operator*=(const double fact) {
  if (active_) bias_ *= fact;
  return *this;
}
// -----------------------------------------------------------------------------
void ModelBiasCorrection::axpy(const double fact, const ModelBiasCorrection & rhs) {
  if (active_) bias_ += fact * rhs.bias_;
}
// -----------------------------------------------------------------------------
double ModelBiasCorrection::dot_product_with(const ModelBiasCorrection & rhs) const {
  double zz = 0.0;
  if (active_) zz = bias_ * rhs.bias_;
  return zz;
}
// -----------------------------------------------------------------------------
size_t ModelBiasCorrection::serialSize() const {
  return 1;
}
// -----------------------------------------------------------------------------
void ModelBiasCorrection::serialize(std::vector<double> & vect) const {
  vect.push_back(bias_);
}
// -----------------------------------------------------------------------------
void ModelBiasCorrection::deserialize(const std::vector<double> & vect, size_t & index) {
  bias_ = vect.at(index);
  ++index;
}
// -----------------------------------------------------------------------------
void ModelBiasCorrection::print(std::ostream & os) const {
  if (active_) {os << std::endl << "ModelBiasCorrection = " << bias_;}
}
// -----------------------------------------------------------------------------
}  // namespace lorenz95
