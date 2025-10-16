/*
 * (C) Copyright 2025- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

namespace oops {

template<typename MODEL1, typename MODEL2>
class ModelDataCoupled : public util::Printable {
  typedef GeometryCoupled<MODEL1, MODEL2> GeometryCoupled_;
 public:
  static const oops::Variables defaultVariables() {return oops::Variables();}

  explicit ModelDataCoupled(const GeometryCoupled_ &) {}
  ~ModelDataCoupled() = default;

  const eckit::LocalConfiguration modelData() const {return eckit::LocalConfiguration();}

 private:
  void print(std::ostream & os) const {}
};

// -----------------------------------------------------------------------------

}  // namespace oops
