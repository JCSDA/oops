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

#include "oops/base/VariableChangeBase.h"

#include "oops/qg/QgTraits.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace qg {
  class GeometryQG;
  class StateQG;

// -----------------------------------------------------------------------------
/// QG change of variable

class ChangeVarQG: public oops::VariableChangeBase<QgTraits> {
 public:
  static const std::string classname() {return "qg::ChangeVarQG";}

  ChangeVarQG(const GeometryQG &, const eckit::Configuration &);
  ~ChangeVarQG();

/// Perform transforms
  void changeVar(const StateQG &, StateQG &) const override;
  void changeVarInverse(const StateQG &, StateQG &) const override;

 private:
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_CHANGEVARQG_H_
