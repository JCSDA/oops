/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <array>
#include <vector>

#include "atlas/util/Geometry.h"
#include "atlas/util/Point.h"

namespace stripack {

// Interface for Fortran implementations
extern "C" {
  void stripack_trmesh_f90(const int &, const double[], const double[], const double[],
                           int[], int[], int[], int &, int[], int[], double[]);
  void stripack_trfind_f90(const int &, const double[],
                           const int &, const double[], const double[], const double[],
                           const int[], const int[], const int[],
                           double[], int[]);
};  // extern "C"

class Triangulation {
 public:
  Triangulation(const std::vector<double> & lats, const std::vector<double> & lons);

  /// Returns true if there is a containing triangle; false if point is outside triangulation
  bool containingTriangleAndBarycentricCoords(
      const std::array<double, 3> & coords,
      int guessIndex,
      std::array<int, 3> & vertexIndices,
      std::array<double, 3> & barycentricCoords) const;

 private:
  // Triangulation data to use in C++
  size_t num_nodes_;
  std::vector<double> xs_;
  std::vector<double> ys_;
  std::vector<double> zs_;
  std::vector<size_t> random_permutation_;
  std::vector<size_t> inverse_random_permutation_;

  // Allocations for STRIPACK connectivity data; should not read from C++
  std::vector<int> list_;
  std::vector<int> lptr_;
  std::vector<int> lend_;
  int lnew_;
  // Allocations for STRIPACK workspaces; should not read from C++
  std::vector<int> near_;
  std::vector<int> next_;
  std::vector<double> dist_;
};

}  // namespace stripack
