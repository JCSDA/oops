/*
 * (C) Crown Copyright 2025-2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

namespace oops {

#ifdef NDEBUG
constexpr bool OOPS_BUILD_TYPE_DEBUG = false;
#else
constexpr bool OOPS_BUILD_TYPE_DEBUG = true;
#endif  // ifdef NDEBUG

}  // namespace oops

