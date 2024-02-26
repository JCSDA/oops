/*
 * (C) Copyright 2017-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QG_MODEL_CHANGEVARQG_H_
#define QG_MODEL_CHANGEVARQG_H_

#include <ostream>
#include <string>

#include "oops/qg/GeometryQG.h"
#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace qg {
  class StateQG;

// -----------------------------------------------------------------------------
/// QG change of variable

class ChangeVarQG : public util::Printable {
 public:
  static const std::string classname() {return "qg::ChangeVarQG";}

  ChangeVarQG(const eckit::Configuration &, const GeometryQG &);
  ~ChangeVarQG();

/// Perform transforms
  void changeVar(StateQG &, const oops::Variables &) const;
  void changeVarInverse(StateQG &, const oops::Variables &) const;

 private:
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_CHANGEVARQG_H_
