/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/ObsBias.h"

#include <iostream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "lorenz95/ObsBiasCorrection.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------
ObsBias::ObsBias(const ObsTable &, const eckit::Configuration & conf)
  : bias_(0.0), active_(false), geovars_(std::vector<std::string>{"x"}), hdiags_()
{
  oops::Log::trace() << "ObsBias::ObsBias conf is:" << conf << std::endl;
  if (conf.has("bias")) {
    bias_ = conf.getDouble("bias");
    active_ = true;
    oops::Log::info() << "ObsBias::ObsBias created, bias = " << bias_ << std::endl;
  }
}
// -----------------------------------------------------------------------------
ObsBias::ObsBias(const ObsBias & other, const bool copy)
  : bias_(0.0), active_(other.active_),
    geovars_(other.geovars_), hdiags_(other.hdiags_) {
  if (active_ && copy) bias_ = other.bias_;
}
// -----------------------------------------------------------------------------
ObsBias & ObsBias::operator+=(const ObsBiasCorrection & dx) {
  if (active_) bias_ += dx.value();
  return *this;
}
// -----------------------------------------------------------------------------
ObsBias & ObsBias::operator=(const ObsBias & rhs) {
  if (active_) bias_ = rhs.bias_;
  return *this;
}
// -----------------------------------------------------------------------------
void ObsBias::print(std::ostream & os) const {
  if (active_) {os << "ObsBias = " << bias_;}
}
// -----------------------------------------------------------------------------
}  // namespace lorenz95
