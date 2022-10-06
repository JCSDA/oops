/*
 * (C) Copyright 2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_FIELDSETOPERATIONS_H_
#define OOPS_UTIL_FIELDSETOPERATIONS_H_

#include "atlas/field.h"

namespace util {

// -----------------------------------------------------------------------------

void FieldSetMultiply(atlas::FieldSet & fset, const atlas::FieldSet & mulFset);
void FieldSetDivide(atlas::FieldSet & fset, const  atlas::FieldSet & divFset);
void FieldSetSqrt(atlas::FieldSet & fset);

// -----------------------------------------------------------------------------

}  // namespace util

#endif  // OOPS_UTIL_FIELDSETOPERATIONS_H_
