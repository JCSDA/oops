/*
 * (C) Crown copyright 2021, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <vector>

#include "model/ObsBiasPreconditioner.h"

namespace qg {
  class ObsBias;
// -----------------------------------------------------------------------------

ObsBiasPreconditioner::ObsBiasPreconditioner(const std::array<double, ObsBias::ntypes> & precond) :
    precond_(precond) {}

void ObsBiasPreconditioner::multiply(const ObsBiasIncrement & dx1, ObsBiasIncrement & dx2) const {
    for (unsigned int jj = 0; jj < ObsBias::ntypes; ++jj) {
      if (precond_[jj] > 0.0) {
        dx2[jj] = dx1[jj] * precond_[jj];
      } else {
        dx2[jj] = 0.0;
      }
    }
}

// -----------------------------------------------------------------------------

}  // namespace qg
