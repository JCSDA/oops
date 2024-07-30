/*
 * (C) Copyright 2020 Met Office UK
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_UTIL_ISANYPOINTINVOLUMEINTERIOR_H_
#define TEST_UTIL_ISANYPOINTINVOLUMEINTERIOR_H_

#include <array>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/container/KDTree.h"
#include "eckit/geometry/KPoint.h"
#include "eckit/testing/Test.h"
#include "oops/../test/TestEnvironment.h"
#include "oops/runs/Test.h"
#include "oops/util/Expect.h"
#include "oops/util/IsAnyPointInVolumeInterior.h"
#include "oops/util/sqr.h"

namespace test {

struct KDTreeBuilder {
  struct EmptyPayload {};

  static const int numDims = 3;
  static const int pointGridSize = 3;

  typedef eckit::geometry::KPoint<numDims> KPoint;

  struct TreeTraits {
    typedef KPoint Point;
    typedef EmptyPayload Payload;
  };

  typedef eckit::KDTreeMemory<TreeTraits> KDTree;

  static KDTree buildTree() {
    std::vector<KDTree::Value> values;
    KPoint point;
    for (int i = 0; i < pointGridSize; ++i) {
      point.data()[i] = i;
      for (int j = 0; j < pointGridSize; ++j) {
        point.data()[j] = j;
        for (int k = 0; k < pointGridSize; ++k) {
          point.data()[k] = k;
          values.push_back(KDTree::Value(KDTree::Point(point), EmptyPayload()));
        }
      }
    }

    KDTree tree;
    tree.build(values);
    return tree;
  }
};

CASE("util/IsAnyPointInVolumeInterior/center") {
  KDTreeBuilder::KDTree tree = KDTreeBuilder::buildTree();
  const int numDims = KDTreeBuilder::numDims;
  const int numCylinderBaseDims = numDims - 1;
  const int pointGridSize = KDTreeBuilder::pointGridSize;
  typedef KDTreeBuilder::KPoint KPoint;

  KPoint halfSize(std::array<double, numDims>{1. / 16., 1. / 8., 1 / 4.});

  KPoint point;
  for (int i = 0; i < pointGridSize; ++i) {
    point.data()[i] = i;
    for (int j = 0; j < pointGridSize; ++j) {
      point.data()[j] = j;
      for (int k = 0; k < pointGridSize; ++k) {
        point.data()[k] = k;

        // Lattice points at the centers of volumes
        EXPECT(util::isAnyPointInBoxInterior(tree, point - halfSize, point + halfSize));
        EXPECT(util::isAnyPointInEllipsoidInterior(tree, point - halfSize, point + halfSize));
        EXPECT(util::isAnyPointInCylinderInterior(tree, point - halfSize, point + halfSize,
                                                  numCylinderBaseDims));
      }
    }
  }
}

CASE("util/IsAnyPointInVolumeInterior/shiftedAlongCartesianAxis") {
  KDTreeBuilder::KDTree tree = KDTreeBuilder::buildTree();
  const int numDims = KDTreeBuilder::numDims;
  const int numCylinderBaseDims = numDims - 1;
  const int pointGridSize = KDTreeBuilder::pointGridSize;
  typedef KDTreeBuilder::KPoint KPoint;

  KPoint halfSize(std::array<double, numDims>{1. / 16., 1. / 8., 1 / 4.});

  KPoint point;
  for (int i = 0; i < pointGridSize; ++i) {
    point.data()[i] = i;
    for (int j = 0; j < pointGridSize; ++j) {
      point.data()[j] = j;
      for (int k = 0; k < pointGridSize; ++k) {
        point.data()[k] = k;

        // Volumes shifted along the the +-x, +-y or +-z direction
        // so that the lattice points are just inside their surfaces
        for (int m = 0; m < 2 * numDims; ++m) {
          const int dim = m / 2;
          const int sign = m % 2 ? 1 : -1;
          KPoint center = point;
          center.data()[dim] += sign * 0.99 * halfSize.data()[dim];
          EXPECT(util::isAnyPointInBoxInterior(tree, center - halfSize, center + halfSize));
          EXPECT(util::isAnyPointInEllipsoidInterior(tree, center - halfSize, center + halfSize));
          EXPECT(util::isAnyPointInCylinderInterior(tree, center - halfSize, center + halfSize,
                                                    numCylinderBaseDims));
        }

        // Volumes shifted along the the +-x, +-y or +-z direction
        // so that the lattice points are just outside their surfaces
        for (int m = 0; m < 2 * numDims; ++m) {
          const int dim = m / 2;
          const int sign = m % 2 ? 1 : -1;
          KPoint center = point;
          center.data()[dim] += sign * 1.01 * halfSize.data()[dim];
          EXPECT_NOT(util::isAnyPointInBoxInterior(tree, center - halfSize, center + halfSize));
          EXPECT_NOT(util::isAnyPointInEllipsoidInterior(tree,
                                                         center - halfSize, center + halfSize));
          EXPECT_NOT(util::isAnyPointInCylinderInterior(tree, center - halfSize, center + halfSize,
                                                        numCylinderBaseDims));
        }
      }
    }
  }
}

CASE("util/IsAnyPointInVolumeInterior/shiftedCloseToEllipsoidSurface") {
  KDTreeBuilder::KDTree tree = KDTreeBuilder::buildTree();
  const int numDims = KDTreeBuilder::numDims;
  const int numCylinderBaseDims = numDims - 1;
  const int pointGridSize = KDTreeBuilder::pointGridSize;
  typedef KDTreeBuilder::KPoint KPoint;

  KPoint halfSize(std::array<double, numDims>{1. / 16., 1. / 8., 1 / 4.});

  KPoint point;
  for (int i = 0; i < pointGridSize; ++i) {
    point.data()[i] = i;
    for (int j = 0; j < pointGridSize; ++j) {
      point.data()[j] = j;
      for (int k = 0; k < pointGridSize; ++k) {
        point.data()[k] = k;

        // Volumes shifted so that the lattice points are just inside the ellipsoid surface
        double xyz = 0;
        for (int m = 0; m < numDims; ++m)
          xyz += util::sqr(1.0 / halfSize.x(m));
        xyz = 1.0 / std::sqrt(xyz);
        KPoint offset(std::initializer_list<double>{xyz, xyz, xyz});
        for (int m = 0; m < 1; ++m) {
          offset.data()[0] *= -1;
          for (int n = 0; n < 1; ++n) {
            offset.data()[1] *= -1;
            for (int p = 0; p < 1; ++p) {
              offset.data()[2] *= -1;

              KPoint center = point + offset * 0.99;
              EXPECT(util::isAnyPointInBoxInterior(
                       tree, center - halfSize, center + halfSize));
              EXPECT(util::isAnyPointInEllipsoidInterior(
                       tree, center - halfSize, center + halfSize));
              EXPECT(util::isAnyPointInCylinderInterior(
                       tree, center - halfSize, center + halfSize, numCylinderBaseDims));
            }
          }
        }

        // Volumes shifted so that the lattice points are just outside the ellipsoid surface
        for (int m = 0; m < 1; ++m) {
          offset.data()[0] *= -1;
          for (int n = 0; n < 1; ++n) {
            offset.data()[1] *= -1;
            for (int p = 0; p < 1; ++p) {
              offset.data()[2] *= -1;

              KPoint center = point + offset * 1.01;
              EXPECT(util::isAnyPointInBoxInterior(
                       tree, center - halfSize, center + halfSize));
              EXPECT_NOT(util::isAnyPointInEllipsoidInterior(
                           tree, center - halfSize, center + halfSize));
              EXPECT(util::isAnyPointInCylinderInterior(
                       tree, center - halfSize, center + halfSize, numCylinderBaseDims));
            }
          }
        }
      }
    }
  }
}

CASE("util/IsAnyPointInVolumeInterior/shiftedCloseToCylinderSurface") {
  KDTreeBuilder::KDTree tree = KDTreeBuilder::buildTree();
  const int numDims = KDTreeBuilder::numDims;
  const int numCylinderBaseDims = numDims - 1;
  const int pointGridSize = KDTreeBuilder::pointGridSize;
  typedef KDTreeBuilder::KPoint KPoint;

  KPoint halfSize(std::array<double, numDims>{1. / 16., 1. / 8., 1 / 4.});

  KPoint point;
  for (int i = 0; i < pointGridSize; ++i) {
    point.data()[i] = i;
    for (int j = 0; j < pointGridSize; ++j) {
      point.data()[j] = j;
      for (int k = 0; k < pointGridSize; ++k) {
        point.data()[k] = k;

        // Volumes shifted so that the lattice points are just inside the cylinder surface
        double xyz = 0;
        for (int m = 0; m < numCylinderBaseDims; ++m)
          xyz += util::sqr(1.0 / halfSize.x(m));
        xyz = 1.0 / std::sqrt(xyz);
        KPoint offset(std::initializer_list<double>{xyz, xyz, 0.99 * halfSize.x(2)});
        for (int m = 0; m < 1; ++m) {
          offset.data()[0] *= -1;
          for (int n = 0; n < 1; ++n) {
            offset.data()[1] *= -1;

            KPoint center = point + offset * 0.99;
            EXPECT(util::isAnyPointInBoxInterior(tree, center - halfSize, center + halfSize));
            EXPECT_NOT(util::isAnyPointInEllipsoidInterior(tree,
                                                           center - halfSize, center + halfSize));
            EXPECT(util::isAnyPointInCylinderInterior(tree, center - halfSize, center + halfSize,
                                                      numCylinderBaseDims));
          }
        }

        // Volumes shifted so that the lattice points are just outside the cylinder surface
        for (int m = 0; m < 1; ++m) {
          offset.data()[0] *= -1;
          for (int n = 0; n < 1; ++n) {
            offset.data()[1] *= -1;

            KPoint center = point + offset * 1.01;
            EXPECT(util::isAnyPointInBoxInterior(tree, center - halfSize, center + halfSize));
            EXPECT_NOT(util::isAnyPointInEllipsoidInterior(tree,
                                                           center - halfSize, center + halfSize));
            EXPECT_NOT(util::isAnyPointInCylinderInterior(tree,
                                                          center - halfSize, center + halfSize,
                                                          numCylinderBaseDims));
          }
        }
      }
    }
  }
}

class IsAnyPointInVolumeInterior : public oops::Test {
 public:
  using oops::Test::Test;
 private:
  std::string testid() const override {return "ufo::test::IsAnyPointInVolumeInterior";}

  void register_tests() const override {}
  void clear() const override {}
};

}  // namespace test

#endif  // TEST_UTIL_ISANYPOINTINVOLUMEINTERIOR_H_
