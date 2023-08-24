/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_IDCHANGEVARIABLE_H_
#define LORENZ95_IDCHANGEVARIABLE_H_

#include <ostream>
#include <string>

#include "oops/util/Printable.h"

#include "lorenz95/Resolution.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace lorenz95 {
  class StateL95;

// -------------------------------------------------------------------------------------------------
/// No change of variable

class IdChangeVariable : public util::Printable {
 public:
  static const std::string classname() {return "lorenz95::IdChangeVariable";}

  IdChangeVariable(const eckit::Configuration &, const Resolution &) {}
  void changeVar(StateL95 &, const oops::Variables &) const {}
  void changeVarInverse(StateL95 &, const oops::Variables &) const {}

 private:
  void print(std::ostream & os) const override {os << "IdChangeVariable";}
};

// -------------------------------------------------------------------------------------------------

}  // namespace lorenz95
#endif  // LORENZ95_IDCHANGEVARIABLE_H_
