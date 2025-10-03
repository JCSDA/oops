/*
 * (C) Copyright 2022-2023 UCAR.
 * (C) Crown copyright 2022-2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_HTLMCALCULATOR_H_
#define OOPS_GENERIC_HTLMCALCULATOR_H_

#include <limits>
#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "oops/base/IncrementSet.h"
#include "oops/generic/HtlmEnsemble.h"
#include "oops/generic/HtlmRegularization.h"

namespace oops {

//------------------------------------------------------------------------------

template <typename MODEL>
class HtlmCalculator {
  typedef Geometry<MODEL>                                          Geometry_;
  typedef HtlmEnsemble<MODEL>                                      HtlmEnsemble_;
  typedef Increment<MODEL>                                         Increment_;
  typedef IncrementSet<MODEL>                                      IncrementSet_;

 public:
  HtlmCalculator(const eckit::Configuration &,
                 const Variables &,
                 const Geometry_ &,
                 const atlas::idx_t,
                 const HtlmEnsemble_ &,
                 const std::vector<atlas::idx_t> &);
  void setOfCoeffs(const IncrementSet_ &, const IncrementSet_ &, atlas::FieldSet &) const;

 private:
  const eckit::LocalConfiguration config_;
  const Variables & updateVars_;
  const atlas::idx_t nLocations_;
  const atlas::idx_t nLevels_;
  const atlas::idx_t influenceSize_;
  const atlas::idx_t halfInfluenceSize_;
  const atlas::idx_t ensembleSize_;
  const atlas::idx_t vectorSize_;
  const std::vector<atlas::idx_t> & owned_;
  const atlas::FieldSet rmsVals_;
  mutable Eigen::MatrixXd M_;
  mutable Eigen::VectorXd linearErrorVector_;
  mutable Eigen::BDCSVD<Eigen::MatrixXd> SVD_;
  std::unique_ptr<HtlmRegularization> regularization_;

  void singularValueDecomposition(const atlas::idx_t, const atlas::array::Range &,
                                  const IncrementSet_ &) const;
  void compute(const atlas::idx_t, const atlas::idx_t, const atlas::array::Range &,
               const IncrementSet_ &, atlas::FieldSet &) const;
};

//------------------------------------------------------------------------------

template <typename MODEL>
HtlmCalculator<MODEL>::HtlmCalculator(const eckit::Configuration & config,
                                      const Variables & updateVars,
                                      const Geometry_ & updateGeometry,
                                      const atlas::idx_t influenceSize,
                                      const HtlmEnsemble_ & ensemble,
                                      const std::vector<atlas::idx_t> & owned)
: config_(config), updateVars_(updateVars), nLocations_(updateGeometry.functionSpace().size()),
  nLevels_(updateGeometry.variableSizes(updateVars_)[0]), influenceSize_(influenceSize),
  halfInfluenceSize_(influenceSize_ / 2),
  ensembleSize_(ensemble.size()), vectorSize_(influenceSize_ * updateVars_.size()), owned_(owned),
  rmsVals_(ensemble.getRmsVals(updateVars_, nLevels_)), M_(vectorSize_, ensembleSize_),
  linearErrorVector_(ensembleSize_), SVD_(vectorSize_, vectorSize_, Eigen::ComputeThinU) {
  // Set up regularization, can be empty
  const eckit::LocalConfiguration regConfig = config_.getSubConfiguration("regularization");
  if (!regConfig.has("parts"))  {
    regularization_ = std::make_unique<HtlmRegularization>(regConfig);
  } else {
    Increment_ regularizationIncrement(updateGeometry, updateVars_, util::DateTime());
    atlas::FieldSet regularizationFieldSet = regularizationIncrement.fieldSet().fieldSet();
    regularization_ = std::make_unique<HtlmRegularizationComponentDependent>(
      regConfig, regularizationFieldSet);
  }
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HtlmCalculator<MODEL>::setOfCoeffs(const IncrementSet_ & linearEnsemble,
                                        const IncrementSet_ & linearErrors,
                                        atlas::FieldSet & coeffsFSet) const {
  // Loop over grid points and levels, and at each:
  // - form the preconditioned matrix of influencing components across each ensemble member, M;
  // - compute the singular value decomposition of M(M^T);
  // - for each variable, compute vectors of coefficients and store in coeffsFSet.
  // Generally the regions of influence are centred on the level of interest, expect near the bottom
  // and top, where they are the bottom-most/top-most (influenceSize_) levels. The loop over k is
  // therefore split at the ends to avoid repeating the same SVD multiple times.
  for (auto i : owned_) {
    atlas::array::Range range(0, influenceSize_);
    singularValueDecomposition(i, range, linearEnsemble);
    for (auto k = 0; k < halfInfluenceSize_; k++) {
      compute(i, k, range, linearErrors, coeffsFSet);
    }
    for (auto k = halfInfluenceSize_; k < nLevels_ - halfInfluenceSize_; k++) {
      range = atlas::array::Range(k - halfInfluenceSize_, k + halfInfluenceSize_ + 1);
      singularValueDecomposition(i, range, linearEnsemble);
      compute(i, k, range, linearErrors, coeffsFSet);
    }
    range = atlas::array::Range(nLevels_ - influenceSize_, nLevels_);
    singularValueDecomposition(i, range, linearEnsemble);
    for (auto k = nLevels_ - halfInfluenceSize_; k < nLevels_; k++) {
      compute(i, k, range, linearErrors, coeffsFSet);
    }
  }
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HtlmCalculator<MODEL>::singularValueDecomposition(const atlas::idx_t i,
                                                       const atlas::array::Range & range,
                                                       const IncrementSet_ & linearEnsemble)
const {
  // M is a matrix where each column forms a vector of (no. variables) segments, each of length
  // influenceSize_, and each column is taken from one ensemble member
  for (size_t v = 0; v < updateVars_.size(); v++) {
    auto rms = atlas::array::make_view<double, 1>(rmsVals_[updateVars_[v].name()]).slice(range);
    for (auto m = 0; m < ensembleSize_; m++) {
      const auto values = atlas::array::make_view<double, 2>(
        linearEnsemble[m].fieldSet()[updateVars_[v].name()]).slice(i, range);
      for (auto s = 0; s < influenceSize_; s++) {
        // Values are normalized by typical magnitudes (from rmsVals_) as preconditioning
        M_(v * influenceSize_ + s, m) = values(s) / rms[s];
      }
    }
  }
  SVD_.compute(M_ * M_.transpose());
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HtlmCalculator<MODEL>::compute(const atlas::idx_t i, const atlas::idx_t k,
                                    const atlas::array::Range & range,
                                    const IncrementSet_ & linearErrors,
                                    atlas::FieldSet & coeffsFSet) const {
  for (const auto & var : updateVars_.variables()) {
    // Produce VectorXd of linear errors at var, i, k for each ensemble member
    for (auto m = 0; m < ensembleSize_; m++) {
      linearErrorVector_(m, 0)
        = atlas::array::make_view<double, 2>(linearErrors[m].fieldSet()[var])(i, k);
    }

    // Add regularization value to singular values
    // TODO(Tom): change regularization to use variables
    const Eigen::ArrayXd singVals = SVD_.singularValues().array()
                                    + regularization_->getRegularizationValue(var, i, k);

    // Compute vector of coeffs for var at i, k, assigning directly into FieldSet
    const double tol = (singVals.matrix().norm() * std::numeric_limits<double>::epsilon());
    Eigen::Map<Eigen::VectorXd> coeffsMap(
      &atlas::array::make_view<double, 3>(coeffsFSet[var])(i, k, 0), vectorSize_);
    coeffsMap.noalias() = SVD_.matrixU()
      * ((singVals > tol).select(singVals.inverse(), 0.0)).matrix().asDiagonal()
      * SVD_.matrixU().transpose() * M_ * linearErrorVector_;
    // (equation 26 in https://doi.org/10.1175/MWR-D-20-0088.1)

    // Un-normalize to account for preconditioning
    for (size_t v = 0; v < updateVars_.size(); v++) {
      auto rms = atlas::array::make_view<double, 1>(rmsVals_[updateVars_[v].name()]).slice(range);
      for (auto s = 0; s < influenceSize_; s++) {
        coeffsMap(v * influenceSize_ + s, 0) /= rms[s];
      }
    }
  }
}

}  // namespace oops

#endif  // OOPS_GENERIC_HTLMCALCULATOR_H_
