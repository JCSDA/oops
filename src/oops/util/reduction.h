/*
 * (C) Copyright 2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#pragma once

namespace atlas {
class Field;
}  // namespace atlas

namespace util {
enum class ExecutionPattern;
}  // namespace util

namespace util {

// Dot product of two fields, local to the calling task.
// This overload gives the caller control over the execution pattern.
double dot_product_on_task(const ExecutionPattern pattern,
    const atlas::Field & field1, const atlas::Field & field2);

// Dot product of two fields, local to the calling task.
// This overload uses the default execution pattern.
double dot_product_on_task(const atlas::Field & field1, const atlas::Field & field2);

}  // namespace util
