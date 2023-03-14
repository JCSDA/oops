/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "oops/coupled/AuxCoupledModel.h"
#include "oops/coupled/GeometryCoupled.h"
#include "oops/coupled/StateCoupled.h"
#include "oops/coupled/VariableChangeCoupled.h"

namespace oops {

template <typename MODEL1, typename MODEL2>
struct TraitCoupled {
 public:
  static std::string name() {return "Coupled";}

  typedef GeometryCoupled<MODEL1, MODEL2>           Geometry;
  // Increment is not used anywhere, but required to build the code
  // with GetValues class (H(x) application); for now pointing to State
  typedef StateCoupled<MODEL1, MODEL2>              Increment;
  typedef StateCoupled<MODEL1, MODEL2>              State;
  typedef AuxCoupledModel<MODEL1, MODEL2>           ModelAuxControl;
  typedef VariableChangeCoupled<MODEL1, MODEL2>     VariableChange;
};

}  // namespace oops

