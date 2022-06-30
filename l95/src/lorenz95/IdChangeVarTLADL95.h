/*
 * (C) Copyright 2021 UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_IDCHANGEVARTLADL95_H_
#define LORENZ95_IDCHANGEVARTLADL95_H_

#include <ostream>
#include <string>

#include "oops/util/Printable.h"
#include "IdChangeVarTLADL95Params.h"
#include "L95Traits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace lorenz95 {
  class Resolution;
  class StateL95;
  class IncrementL95;

// -----------------------------------------------------------------------------
/// L95 linear change of variable

class IdChangeVarTLADL95: public util::Printable {
 public:
  typedef IdChangeVarTLADL95Params Parameters_;
  static const std::string classname() {return "lorenz95::IdChangeVarTLADL95";}

  explicit IdChangeVarTLADL95(const Resolution &, const Parameters_ &) {}
  ~IdChangeVarTLADL95() {}

/// Perform linear transforms
  void multiply(IncrementL95 &, const oops::Variables &) const {}
  void multiplyInverse(IncrementL95 &, const oops::Variables &) const {}
  void multiplyAD(IncrementL95 &, const oops::Variables &) const {}
  void multiplyInverseAD(IncrementL95 &, const oops::Variables &) const {}

  void setTrajectory(const StateL95 &, const StateL95 &) {}

 private:
  void print(std::ostream & os) const {os << "IdChangeVarTLADL95";}
};
// -----------------------------------------------------------------------------

}  // namespace lorenz95
#endif  // LORENZ95_IDCHANGEVARTLADL95_H_
