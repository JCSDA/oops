/*
 * (C) Copyright 2020 UCAR.
 * (C) Crown Copyright 2024, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/assimilation/MinimizerUtils.h"

#include "oops/util/formats.h"
#include "oops/util/Logger.h"

namespace oops {

void printNormReduction(int iteration, const double & grad, const double & norm) {
  Log::info() << "  Residual norm (" << std::setw(2) << iteration << ") = "
              << util::full_precision(grad) << std::endl
              << "  Norm reduction (" << std::setw(2) << iteration << ") = "
              << util::full_precision(norm) << std::endl << std::endl;
}

void printQuadraticCostFunction(int iteration, const double & costJ,
                                const double & costJb, const double & costJoJc) {
  Log::info() << "  Quadratic cost function: J   (" << std::setw(2) << iteration << ") = "
              << util::full_precision(costJ)        << std::endl
              << "  Quadratic cost function: Jb  (" << std::setw(2) << iteration << ") = "
              << util::full_precision(costJb)       << std::endl
              << "  Quadratic cost function: JoJc(" << std::setw(2) << iteration << ") = "
              << util::full_precision(costJoJc)     << std::endl << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops
