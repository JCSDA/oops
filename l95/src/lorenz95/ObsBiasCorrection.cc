/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/ObsBiasCorrection.h"

#include <iostream>
#include <string>

#include "util/Logger.h"
#include "lorenz95/ObsBias.h"
#include "lorenz95/ObsBiasCovariance.h"
#include "eckit/config/Configuration.h"


using oops::Log;

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------
ObsBiasCorrection::ObsBiasCorrection(const eckit::Configuration & conf)
  : bias_(0.0), active_(false)
{
  active_ = conf.has("standard_deviation");
  if (active_) {Log::trace() << "ObsBiasCorrection::ObsBiasCorrection created." << std::endl;}
}
// -----------------------------------------------------------------------------
ObsBiasCorrection::ObsBiasCorrection(const ObsBiasCorrection & other,
                                     const bool copy)
  : bias_(0.0), active_(other.active_)
{
  if (active_ && copy) bias_ = other.bias_;
}
// -----------------------------------------------------------------------------
ObsBiasCorrection::ObsBiasCorrection(const ObsBiasCorrection & other,
                                     const eckit::Configuration &)
  : bias_(0.0), active_(other.active_)
{
  if (active_) bias_ = other.bias_;
}
// -----------------------------------------------------------------------------
void ObsBiasCorrection::diff(const ObsBias & b1, const ObsBias & b2) {
  if (active_) bias_ = b1.value() - b2.value();
}
// -----------------------------------------------------------------------------
void ObsBiasCorrection::zero() {
  bias_ = 0.0;
}
// -----------------------------------------------------------------------------
ObsBiasCorrection & ObsBiasCorrection::operator=(const ObsBiasCorrection & rhs) {
  if (active_) {
    bias_ = rhs.bias_;
  } else {
    bias_ = 0.0;
  }
  return *this;
}
// -----------------------------------------------------------------------------
ObsBiasCorrection & ObsBiasCorrection::operator+=(const ObsBiasCorrection & rhs) {
  if (active_) bias_ += rhs.bias_;
  return *this;
}
// -----------------------------------------------------------------------------
ObsBiasCorrection & ObsBiasCorrection::operator-=(const ObsBiasCorrection & rhs) {
  if (active_) bias_ -= rhs.bias_;
  return *this;
}
// -----------------------------------------------------------------------------
ObsBiasCorrection & ObsBiasCorrection::operator*=(const double fact) {
  if (active_) bias_ *= fact;
  return *this;
}
// -----------------------------------------------------------------------------
void ObsBiasCorrection::axpy(const double fact, const ObsBiasCorrection & rhs) {
  if (active_) bias_ += fact * rhs.bias_;
}
// -----------------------------------------------------------------------------
double ObsBiasCorrection::dot_product_with(const ObsBiasCorrection & rhs) const {
  double zz = 0.0;
  if (active_) zz = bias_ * rhs.bias_;
  return zz;
}
// -----------------------------------------------------------------------------
void ObsBiasCorrection::print(std::ostream & os) const {
  if (active_) {os << std::endl << "ObsBiasCorrection = " << bias_;}
}
// -----------------------------------------------------------------------------
}  // namespace lorenz95
