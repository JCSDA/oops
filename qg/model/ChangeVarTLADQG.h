/*
 * (C) Copyright 2017-2018  UCAR.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QG_MODEL_CHANGEVARTLADQG_H_
#define QG_MODEL_CHANGEVARTLADQG_H_

#include <ostream>
#include <string>

#include "oops/util/Printable.h"

#include "oops/qg/QgFortran.h"

// Forward declarations
namespace eckit {
  class Configuration;
}

namespace qg {
  class GeometryQG;
  class StateQG;
  class IncrementQG;

// -----------------------------------------------------------------------------
/// QG linear change of variable

class ChangeVarTLADQG: public util::Printable {
 public:
  static const std::string classname() {return "qg::ChangeVarQG";}

  ChangeVarTLADQG(const StateQG &, const StateQG &, const GeometryQG &,
                  const eckit::Configuration &);
  ~ChangeVarTLADQG();

/// Perform linear transforms
  void multiply(const IncrementQG &, IncrementQG &) const;
  void multiplyInverse(const IncrementQG &, IncrementQG &) const;
  void multiplyAD(const IncrementQG &, IncrementQG &) const;
  void multiplyInverseAD(const IncrementQG &, IncrementQG &) const;

 private:
  void print(std::ostream &) const override;

// Data
  F90chvar keyConfig_;
};
// -----------------------------------------------------------------------------

}  // namespace qg
#endif  // QG_MODEL_CHANGEVARTLADQG_H_
