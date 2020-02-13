/*
 * (C) Copyright 2020 UK Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_UTIL_ISANYPOINTINVOLUMEINTERIOR_H_
#define OOPS_UTIL_ISANYPOINTINVOLUMEINTERIOR_H_

#include "eckit/container/KDTree.h"

#include "oops/util/IsPointInVolumeInterior.h"

namespace util {

/// \brief Class used to test whether any point in a kd-tree lies in the interior of a particular
/// volume.
template <typename TreeTraits>
class IsAnyPointInVolumeInterior {
 public:
  typedef eckit::KDTreeX<TreeTraits> KDTree;
  typedef typename KDTree::Point Point;

 private:
  typedef typename KDTree::Alloc Alloc;
  typedef typename KDTree::Node Node;

 public:
  /// \brief Returns true if any point in \p tree lies in the interior of a specified volume.
  ///
  /// \param tree
  ///   Tree to search.
  /// \param bboxLbound
  ///   Lower-left corner of the axis-aligned bounding box of the volume.
  /// \param bboxUbound
  ///   Upper-right corner of the axis-aligned bounding box of the volume.
  /// \param isPointInVolumeInterior
  ///   A functor taking a point as parameter and returning true if that point lies in the interior
  ///   of the volume of interest and false otherwise.
  template <typename Predicate>
  static bool isAnyPointInVolumeInterior(const KDTree &tree,
                                         const Point &bboxLbound, const Point &bboxUbound,
                                         const Predicate &isPointInVolumeInterior) {
    if (!tree.root_) {
      return false;
    }
    Alloc &alloc = tree.alloc_;
    Node *root = alloc.convert(tree.root_, static_cast<Node*>(nullptr));

    return isAnyPointInVolumeInterior(root, alloc, bboxLbound, bboxUbound,
                                      isPointInVolumeInterior);
  }

 private:
  template <typename Predicate>
  static bool isAnyPointInVolumeInterior(const Node *node, Alloc &alloc,
                                         const Point &bboxLbound, const Point &bboxUbound,
                                         const Predicate &isPointInVolumeInterior) {
    const Point &point = node->value().point();

    if (isPointInVolumeInterior(point))
      return true;

    const size_t axis = node->axis();

    if (bboxLbound.x(axis) < point.x(axis))
      if (Node *left = node->left(alloc))
        if (isAnyPointInVolumeInterior(left, alloc, bboxLbound, bboxUbound,
                                       isPointInVolumeInterior))
          return true;

    if (bboxUbound.x(axis) > point.x(axis))
      if (Node *right = node->right(alloc))
        if (isAnyPointInVolumeInterior(right, alloc, bboxLbound, bboxUbound,
                                       isPointInVolumeInterior))
          return true;

    return false;
  }
};

/// \brief Returns true if any point in \p tree is in the interior of the axis-aligned box
/// with bottom-left and top-right corners at \p boxLbound and \p boxUbound.
template <typename TreeTraits>
bool isAnyPointInBoxInterior(const eckit::KDTreeX<TreeTraits> &tree,
                             const typename eckit::KDTreeX<TreeTraits>::Point &boxLbound,
                             const typename eckit::KDTreeX<TreeTraits>::Point &boxUbound) {
  typedef typename eckit::KDTreeX<TreeTraits>::Point Point;
  return IsAnyPointInVolumeInterior<TreeTraits>::isAnyPointInVolumeInterior(
        tree, boxLbound, boxUbound,
        [&boxLbound, &boxUbound](const Point &point) {
          return isPointInBoxInterior(point, boxLbound, boxUbound);
  });
}

/// \brief Returns true if any point in \p tree is in the interior of the axis-aligned ellipsoid
/// whose minimum axis-aligned bounding box is the box with bottom-left and top-right corners
/// \p bboxLbound and \p bboxUbound.
template <typename TreeTraits>
bool isAnyPointInEllipsoidInterior(
    const eckit::KDTreeX<TreeTraits> &tree,
    const typename eckit::KDTreeX<TreeTraits>::Point &bboxLbound,
    const typename eckit::KDTreeX<TreeTraits>::Point &bboxUbound) {
  typedef typename eckit::KDTreeX<TreeTraits>::Point Point;
  return IsAnyPointInVolumeInterior<TreeTraits>::isAnyPointInVolumeInterior(
        tree, bboxLbound, bboxUbound,
        [&bboxLbound, &bboxUbound](const Point &point) {
          return isPointInEllipsoidInterior(point, bboxLbound, bboxUbound);
  });
}

/// \brief Returns true if any point in \p tree is in the interior of a cylinder or its
/// n-dimensional generalization.
///
/// Let
///   * n := eckit::KDTreeX<TreeTraits>::Point::DIMS
///   * m := numCylinderBaseDims
///   * C := (bboxLbound + bboxUbound)/2
///   * R := (bboxUbound - bboxLbound)/2.
///
/// Then the cylinder in question is defined as the Cartesian product of the axis-oriented
/// m-dimensional ellipsoid with center at C[:m] and semi-axes of length R[:m] and the
/// axis-oriented (n-m)-dimensional box with bottom-left corner at bboxLbound[m:] and top-right
/// corner at bboxUbound[m:].
template <typename TreeTraits>
bool isAnyPointInCylinderInterior(const eckit::KDTreeX<TreeTraits> &tree,
                                  const typename eckit::KDTreeX<TreeTraits>::Point &bboxLbound,
                                  const typename eckit::KDTreeX<TreeTraits>::Point &bboxUbound,
                                  size_t numCylinderBaseDims) {
  typedef typename eckit::KDTreeX<TreeTraits>::Point Point;
  return IsAnyPointInVolumeInterior<TreeTraits>::isAnyPointInVolumeInterior(
        tree, bboxLbound, bboxUbound,
        [&bboxLbound, &bboxUbound, numCylinderBaseDims](const Point &point) {
          return isPointInCylinderInterior(point, bboxLbound, bboxUbound, numCylinderBaseDims);
  });
}
}  // namespace util

#endif  // OOPS_UTIL_ISANYPOINTINVOLUMEINTERIOR_H_

