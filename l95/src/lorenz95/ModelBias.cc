/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Crown Copyright 2023, the Met Office.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/ModelBias.h"

#include <iostream>
#include <string>

#include "eckit/config/Configuration.h"
#include "lorenz95/ModelBiasCorrection.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------
ModelBias::ModelBias(const Resolution &, const eckit::Configuration & conf)
  : bias_(0.0), active_(false)
{
  oops::Log::trace() << "ModelBias::ModelBias conf is:" << conf << std::endl;
  if (conf.has("bias")) {
    bias_ = conf.getDouble("bias");
    active_ = true;
  }
}
// -----------------------------------------------------------------------------
ModelBias::ModelBias(const Resolution &, const ModelBias & other)
  : bias_(0.0), active_(other.active_)
{
  if (active_) bias_ = other.bias_;
}
// -----------------------------------------------------------------------------
ModelBias::ModelBias(const ModelBias & other, const bool copy)
  : bias_(0.0), active_(other.active_)
{
  if (active_ && copy) bias_ = other.bias_;
}
// -----------------------------------------------------------------------------
ModelBias & ModelBias::operator=(const ModelBias & other) {
  if (active_) bias_ = other.bias_;
  return *this;
}
// -----------------------------------------------------------------------------
ModelBias & ModelBias::operator+=(const ModelBiasCorrection & dx) {
  if (active_) bias_ += dx.bias();
  return *this;
}
// -----------------------------------------------------------------------------
void ModelBias::print(std::ostream & os) const {
  if (active_) {os << "ModelBias = " << bias_;}
}
// -----------------------------------------------------------------------------
}  // namespace lorenz95
