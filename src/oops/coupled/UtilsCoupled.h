/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <vector>

namespace oops {

class Variables;

/// Returns \p varsToSplit variables split by coupled model components according to which
/// variables are available for different model components.
/// The return vector has the same size as \p varsOfCoupledModels. Each of the element of
/// the return vector should only contain variables that also exist in the corresponding
/// element of \p varsOfCoupledModels. The superset of all elements in the return vector
/// should be the same as \p varsToSplit. If not all variables in \p varsToSplit are
/// available in model components, an exception is thrown.
std::vector<Variables> splitVariables(const Variables & varsToSplit,
                                      const std::vector<Variables> & varsOfCoupledModels);

}  // namespace oops
