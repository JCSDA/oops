/*
 * (C) Copyright 2024 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>

#include "oops/base/Variables.h"
#include "oops/util/Printable.h"

namespace lorenz95 {
  class Resolution;

// -------------------------------------------------------------------------------------------------

class ModelData : public util::Printable {
 public:
  static const std::string classname() {return "lorenz95::ModelData";}
  static const oops::Variables defaultVariables();

  explicit ModelData(const Resolution &);
  ~ModelData();

  const eckit::LocalConfiguration modelData() const;

 private:
  void print(std::ostream &) const override;
};

// -------------------------------------------------------------------------------------------------

}  // namespace lorenz95
