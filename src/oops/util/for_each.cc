/*
 * (C) Copyright 2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 *
 */

#include "oops/util/for_each.h"

namespace util {
namespace details {

ExecutionPattern getDefaultForEachExecutionPattern() { return ExecutionPattern::parallel; }

IndexRange getDefaultForEachIndexRange() { return IndexRange::exclude_halo; }

}  // namespace details
}  // namespace util
