/*
 * (C) Copyright 2009-2016 ECMWF.
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
ModelBias::ModelBias(const Resolution &, const Parameters_ & parameters)
  : bias_(0.0), active_(false)
{
  oops::Log::trace() << "ModelBias::ModelBias conf is:" << parameters << std::endl;
  if (parameters.bias.value() != boost::none) {
    bias_ = *parameters.bias.value();
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
