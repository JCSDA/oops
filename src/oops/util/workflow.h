/*
 * (C) Copyright 2021-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

namespace util {

// -----------------------------------------------------------------------------

void use_ecflow();
void update_workflow_meter(const std::string &, const int);

// -----------------------------------------------------------------------------

}  // namespace util
