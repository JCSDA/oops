/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QG_MODEL_OBSBIASPRECONDITIONER_H_
#define QG_MODEL_OBSBIASPRECONDITIONER_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "model/ObsBias.h"
#include "model/ObsBiasIncrement.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace qg {
// -----------------------------------------------------------------------------

/// Class for VarBC preconditioner.

class ObsBiasPreconditioner : public util::Printable,
                              private boost::noncopyable,
                              private util::ObjectCounter<ObsBiasPreconditioner> {
 public:
  static const std::string classname() {return "qg::ObsBiasPreconditioner";}

/// Constructor, destructor
  explicit ObsBiasPreconditioner(const std::array<double, ObsBias::ntypes> &);

/// Linear algebra operators
  void multiply(const ObsBiasIncrement &, ObsBiasIncrement &) const;

 private:
  void print(std::ostream &) const {}
  const std::array<double, ObsBias::ntypes> precond_;
};

// -----------------------------------------------------------------------------

}  // namespace qg

#endif  // QG_MODEL_OBSBIASPRECONDITIONER_H_
