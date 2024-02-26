/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_OBSBIASPRECONDITIONER_H_
#define LORENZ95_OBSBIASPRECONDITIONER_H_

#include <ostream>
#include <string>
#include <boost/noncopyable.hpp>

#include "lorenz95/ObsBiasCorrection.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace lorenz95 {
// -----------------------------------------------------------------------------

/// Class for VarBC preconditioner.

class ObsBiasPreconditioner : public util::Printable,
                          private boost::noncopyable,
                          private util::ObjectCounter<ObsBiasPreconditioner> {
 public:
  static const std::string classname() {return "lorenz95::ObsBiasPreconditioner";}

/// Constructor
  explicit ObsBiasPreconditioner(const double &);

/// Linear algebra operators
  void multiply(const ObsBiasCorrection & dx1, ObsBiasCorrection & dx2) const;

 private:
  void print(std::ostream &) const {}
  const double precond_;
};

// -----------------------------------------------------------------------------

}  // namespace lorenz95

#endif  // LORENZ95_OBSBIASPRECONDITIONER_H_
