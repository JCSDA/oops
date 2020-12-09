/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_MINIMIZERUTILS_H_
#define OOPS_ASSIMILATION_MINIMIZERUTILS_H_

namespace oops {

/// Prints to Log::info gradient reduction \p grad and normalized gradient reduction \p norm
/// for iteration \p iteration
void printNormReduction(int iteration, const double & grad, const double & norm);
/// Prints to Log::info cost function values for \p costJ, \p costJb, \p costJoJc for
/// iteration \p iteration
void printQuadraticCostFunction(int iteration, const double & costJ, const double & costJb,
                                const double & costJoJc);
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_MINIMIZERUTILS_H_
