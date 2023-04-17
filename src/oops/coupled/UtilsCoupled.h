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
/// \p varsOfCoupledModels specifies which variables are available for which model
/// components. A variable name can only be in one of the elements of \p varsOfCoupledModels.
/// If a variable name is in more than one elements of \p varsOfCoupledModels, an exception
/// is thrown.
/// The return vector has the same size as \p varsOfCoupledModels. Each of the element of
/// the return vector should only contain variables that also exist in the corresponding
/// element of \p varsOfCoupledModels. The superset of all elements in the return vector
/// should be the same as \p varsToSplit. If not all variables in \p varsToSplit are
/// available in model components, an exception is thrown.
std::vector<Variables> splitVariables(const Variables & varsToSplit,
                                      const std::vector<Variables> & varsOfCoupledModels);

}  // namespace oops
