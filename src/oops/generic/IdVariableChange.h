/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_IDVARIABLECHANGE_H_
#define OOPS_GENERIC_IDVARIABLECHANGE_H_

#include <ostream>
#include <string>

#include "oops/base/VariableChangeBase.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/State.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// No change of variable

template <typename MODEL>
class IdVariableChange : public VariableChangeBase<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;
 public:
  static const std::string classname() {return "oops::IdVariableChange";}

  IdVariableChange(const Geometry_ &, const eckit::Configuration & conf)
    : VariableChangeBase<MODEL>(conf) {}
  virtual ~IdVariableChange() {}

/// Perform identity change of variable
  void changeVar(const State_ & x1, State_ & x2) const override {x2 = x1;}
  void changeVarInverse(const State_ & x1, State_ & x2) const override {x2 = x1;}

 private:
  void print(std::ostream &) const override {}
};

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_GENERIC_IDVARIABLECHANGE_H_
