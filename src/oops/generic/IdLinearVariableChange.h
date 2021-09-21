/*
 * (C) Copyright 2017-2018  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_IDLINEARVARIABLECHANGE_H_
#define OOPS_GENERIC_IDLINEARVARIABLECHANGE_H_

#include <ostream>
#include <string>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/LinearVariableChangeBase.h"
#include "oops/base/State.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {

// -----------------------------------------------------------------------------
/// No change of variable

template <typename MODEL>
class IdLinearVariableChange : public LinearVariableChangeBase<MODEL> {
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef State<MODEL>               State_;
 public:
  static const std::string classname() {return "oops::IdLinearVariableChange";}

  IdLinearVariableChange(const State_ &, const State_ &, const Geometry_ &,
                  const eckit::Configuration & conf): LinearVariableChangeBase<MODEL>(conf) {}
  virtual ~IdLinearVariableChange() {}

/// Perform linear transforms
  void multiply(const Increment_ & dx1, Increment_ & dx2) const override {dx2 = dx1;}
  void multiplyInverse(const Increment_ & dx1, Increment_ & dx2) const override {dx2 = dx1;}
  void multiplyAD(const Increment_ & dx1, Increment_ & dx2) const override {dx2 = dx1;}
  void multiplyInverseAD(const Increment_ & dx1, Increment_ & dx2) const override {dx2 = dx1;}

 private:
  void print(std::ostream & os) const override {os << "IdVariableChange";}
};

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_GENERIC_IDLINEARVARIABLECHANGE_H_
