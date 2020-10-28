/*
 * (C) Copyright 2020 UK Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_ISPOINTINVOLUMEINTERIOR_H_
#define OOPS_UTIL_ISPOINTINVOLUMEINTERIOR_H_

#include "eckit/container/KDTree.h"
#include "oops/util/sqr.h"

namespace util {

/// \brief Returns true if \p point is in the interior of the axis-aligned box
/// with bottom-left and top-right corners at \p lbound and \p ubound.
template <typename Point>
bool isPointInBoxInterior(
    const Point &point, const Point &lbound, const Point &ubound) {
  for (size_t d = 0; d < Point::DIMS; ++d) {
    if (point.x(d) <= lbound.x(d) || point.x(d) >= ubound.x(d))
      return false;
  }
  return true;
}

/// \brief Returns true if \p point is in the interior of the axis-aligned ellipsoid
/// whose minimum axis-aligned bounding box is the box with bottom-left and top-right corners
/// \p lbound and \p ubound.
template <typename Point>
bool isPointInEllipsoidInterior(
    const Point &point, const Point &lbound, const Point &ubound) {
  if (!isPointInBoxInterior(point, lbound, ubound))
    return false;

  double lhs = 0;
  for (size_t d = 0; d < Point::DIMS; ++d) {
    double centre = 0.5 * (ubound.x(d) + lbound.x(d));
    double radius = 0.5 * (ubound.x(d) - lbound.x(d));
    lhs += sqr((point.x(d) - centre) / radius);
    if (lhs >= 1.0)
      return false;
  }

  return true;
}

/// \brief Returns true if \p point is in the interior of a cylinder or its n-dimensional
/// generalization.
///
/// Let
///   * n := point::DIMS
///   * m := numCylinderBaseDims
///   * C := (lbound + ubound)/2
///   * R := (ubound - lbound)/2.
///
/// Then the cylinder in question is defined as the Cartesian product of the axis-oriented
/// m-dimensional ellipsoid with center at C[:m] and semi-axes of length R[:m] and the
/// (n-m)-dimensional axis-oriented box with bottom-left corner at lbound[m:] and top-right corner
/// at ubound[m:].
template <typename Point>
bool isPointInCylinderInterior(
    const Point &point, const Point &lbound, const Point &ubound, size_t numCylinderBaseDims) {
  if (!isPointInBoxInterior(point, lbound, ubound))
    return false;

  double lhs = 0;
  for (size_t d = 0; d < numCylinderBaseDims; ++d) {
      double centre = 0.5 * (ubound.x(d) + lbound.x(d));
      double radius = 0.5 * (ubound.x(d) - lbound.x(d));
      lhs += sqr((point.x(d) - centre) / radius);
      if (lhs >= 1.0)
        return false;
  }

  return true;
}

}  // namespace util

#endif  // OOPS_UTIL_ISPOINTINVOLUMEINTERIOR_H_

