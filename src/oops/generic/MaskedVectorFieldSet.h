// (C) Crown Copyright 2026 Met Office
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#pragma once

#include <vector>

#include "atlas/field/FieldSet.h"

namespace oops {

/// \brief Read \p fieldSet's values into \p dataVector, restricted to
/// \p locationMask.
///
/// \p dataVector is resized to fit \p fieldSet and \p locationMask. The
/// vector is laid out as [field][masked location][level], following
/// \p fieldSet's own field order and level counts.
void fieldSetToVector(const std::vector<bool>& locationMask,
                     const atlas::FieldSet& fieldSet,
                     std::vector<double>& dataVector);

/// \brief Write \p dataVector's values into \p fieldSet, restricted to
/// \p locationMask.
///
/// \p dataVector's size must already match \p fieldSet and \p locationMask,
/// or eckit::BadValue is thrown. The vector is laid out as [field][masked
/// location][level], following \p fieldSet's own field order and level
/// counts.
void vectorToFieldSet(const std::vector<bool>& locationMask,
                     const std::vector<double>& dataVector,
                     atlas::FieldSet& fieldSet);

}  // namespace oops

