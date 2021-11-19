/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "lorenz95/ObsBiasPreconditioner.h"

namespace lorenz95 {

// -----------------------------------------------------------------------------

ObsBiasPreconditioner::ObsBiasPreconditioner(const double & precond) :
    precond_(precond) {}

void ObsBiasPreconditioner::multiply(const ObsBiasCorrection & dx1, ObsBiasCorrection & dx2) const {
  dx2 = dx1;
  dx2 *= precond_;
}



// -----------------------------------------------------------------------------

}  // namespace lorenz95
