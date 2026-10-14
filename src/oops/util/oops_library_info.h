/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#pragma once

#include <ostream>

namespace util {
void oops_library_info(std::ostream & os);
}  // namespace util

extern "C" {
const char * oops_version();
const char * oops_git_sha1();
}

