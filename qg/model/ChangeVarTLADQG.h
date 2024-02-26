/*
 * (C) Copyright 2017-2021  UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QG_MODEL_CHANGEVARTLADQG_H_
#define QG_MODEL_CHANGEVARTLADQG_H_

#include <ostream>
#include <string>

#include "oops/util/Printable.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace qg {
  class GeometryQG;
  class StateQG;
  class IncrementQG;

// -----------------------------------------------------------------------------
/// QG linear change of variable

class ChangeVarTLADQG: public util::Printable {
 public:
  static const std::string classname() {return "qg::ChangeVarTLADQG";}

  ChangeVarTLADQG(const GeometryQG &, const eckit::Configuration &);
  ~ChangeVarTLADQG();

/// Perform linear transforms
  void changeVarTL(IncrementQG &, const oops::Variables &) const;
  void changeVarInverseTL(IncrementQG &, const oops::Variables &) const;
  void changeVarAD(IncrementQG &, const oops::Variables &) const;
  void changeVarInverseAD(IncrementQG &, const oops::Variables &) const;

  void changeVarTraj(const StateQG &, const oops::Variables &);

 private:
  void print(std::ostream &) const override;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_CHANGEVARTLADQG_H_
