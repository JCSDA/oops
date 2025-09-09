/*
 * (C) Crown Copyright 2025, the Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_SUBENSEMBLESPLITTER_H_
#define OOPS_ASSIMILATION_SUBENSEMBLESPLITTER_H_

#include <Eigen/Sparse>
#include <string>
#include <tuple>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

namespace oops {

/// Subensemble splitter class
/*!
 * Class that partitions ensemble members into subensembles.
 *
 * It can generate projection matrices that can project
 * out members of subensembles needed for cross validation.
 */

// -----------------------------------------------------------------------------

class SubensembleSplitter {
 public:
  /// \brief Create unmodulated subensemble splitter for ensemble
  /// \param[in]  nens     No. ensemble members
  /// \param[in]  config   Configuration options
  SubensembleSplitter(const size_t, const eckit::LocalConfiguration &);

  /// \brief Create subensemble splitter for ensemble
  /// \param[in]  nens     No. ensemble members
  /// \param[in]  neig     No. eigenvalues per member
  /// \param[in]  config   Configuration options
  SubensembleSplitter(const size_t, const size_t, const eckit::LocalConfiguration &);

  ~SubensembleSplitter() {}

  /// \brief Splitting method based on configuration option
  void split();

  /// \brief Gets projection matrices and excluded members for a given excluded subensemble
  /// \param[in]   xclsub              Subensemble to exclude
  /// \param[in]   modulated           Boolean to trigger whether or not to get the
  ///                                  excluded projection matrix for modulated ensemble
  /// \param[out]  excludedProjection  Projection matrix that projects out the excluded members
  /// \param[out]  includedProjection  Projection matrix that projects out the included members
  /// \param[out]  excludedMembers     Vector of excluded members, used for testing
  const std::tuple<Eigen::SparseMatrix<float>,
                   Eigen::SparseMatrix<float>,
                   std::vector<size_t>> getProjectionMatrices(const size_t, const bool);

 private:
  // Methods
  void ctorHelper(const eckit::LocalConfiguration &);
  void contiguousPartition(const std::vector<size_t> &);

  // Data
  size_t nens_;
  size_t nsubens_;
  size_t neig_;
  std::string splittingMethod_;
  size_t seed_;
  size_t subsize_;
  std::vector<size_t> subens_;
};

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_SUBENSEMBLESPLITTER_H_
