/*
 * (C) Copyright 2022-2023 UCAR.
 * (C) Crown copyright 2022-2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_HTLMCALCULATOR_H_
#define OOPS_GENERIC_HTLMCALCULATOR_H_

#include <memory>
#include <string>
#include <unordered_map>
#include <vector>

#include "oops/base/IncrementEnsemble.h"
#include "oops/generic/HtlmRegularization.h"

namespace oops {

template <typename MODEL>
class HtlmCalculatorParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(HtlmCalculatorParameters, Parameters);

 public:
  Parameter<HtlmRegularizationParameters> regularization{"regularization", {}, this};
  Parameter<bool> rmsScaling{"rms scaling", true, this};
};

//------------------------------------------------------------------------------

template <typename MODEL>
class HtlmCalculator {
  typedef Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>    EigenMatrix;
  typedef Eigen::Matrix<double, Eigen::Dynamic, 1>                 EigenVector;
  typedef Geometry<MODEL>                                          Geometry_;
  typedef HtlmCalculatorParameters<MODEL>                          Parameters_;
  typedef Increment<MODEL>                                         Increment_;
  typedef IncrementEnsemble<MODEL>                                 IncrementEnsemble_;

 public:
  HtlmCalculator(const Parameters_ &,
                 const Variables &,
                 const Geometry_ &,
                 const atlas::idx_t,
                 const atlas::idx_t);
  void setOfCoeffs(const IncrementEnsemble_ &, const IncrementEnsemble_ &, atlas::FieldSet &) const;

 private:
  const Parameters_ & params_;
  const Variables & updateVars_;
  const atlas::idx_t nLocations_;
  const atlas::idx_t nLevels_;
  const atlas::idx_t influenceSize_;
  const atlas::idx_t halfInfluenceSize_;
  const atlas::idx_t nLevelsMinusInfluenceSize_;
  const atlas::idx_t ensembleSize_;
  const atlas::idx_t vectorSize_;
  mutable std::unordered_map<std::string, std::vector<double>> rmsVals_;
  std::unique_ptr<HtlmRegularization> regularization_;

  void computeVectorsAt(const atlas::idx_t, const atlas::idx_t,
                        const IncrementEnsemble_ &, const IncrementEnsemble_ &,
                        atlas::FieldSet &) const;
  const EigenMatrix makeInfluenceMatrix(const atlas::idx_t, const atlas::idx_t,
                                        const IncrementEnsemble_ &) const;
  const EigenVector computeVector(const Eigen::BDCSVD<EigenMatrix> &,
                                  const EigenMatrix &, const EigenMatrix &, const EigenVector &,
                                  const double &) const;
};

//------------------------------------------------------------------------------

template <typename MODEL>
HtlmCalculator<MODEL>::HtlmCalculator(const Parameters_ & params,
                                      const Variables & updateVars,
                                      const Geometry_ & updateGeometry,
                                      const atlas::idx_t influenceSize,
                                      const atlas::idx_t ensembleSize)
: params_(params), updateVars_(updateVars), nLocations_(updateGeometry.functionSpace().size()),
  nLevels_(updateGeometry.variableSizes(updateVars_)[0]), influenceSize_(influenceSize),
  halfInfluenceSize_(influenceSize_ / 2), nLevelsMinusInfluenceSize_(nLevels_ - influenceSize_),
  ensembleSize_(ensembleSize), vectorSize_(influenceSize_ * updateVars_.size()) {
  // Set up regularization
  if (params_.regularization.value().parts.value() == boost::none) {
    regularization_ = std::make_unique<HtlmRegularization>(params_.regularization.value());
  } else {
    Increment_ regularizationIncrement(updateGeometry, updateVars_, util::DateTime());
    atlas::FieldSet regularizationFieldSet = regularizationIncrement.fieldSet().fieldSet();
    regularization_ = std::make_unique<HtlmRegularizationComponentDependent>(
      params_.regularization.value(), regularizationFieldSet);
  }
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HtlmCalculator<MODEL>::setOfCoeffs(const IncrementEnsemble_ & linearEnsemble,
                                        const IncrementEnsemble_ & linearErrors,
                                        atlas::FieldSet & coeffsFSet) const {
  // Calculate RMS-by-level scaling values if required
  if (params_.rmsScaling) {
    rmsVals_.clear();
    for (const auto & var : updateVars_.variables()) {
      rmsVals_.emplace(var, linearEnsemble[0].rmsByLevel(var));
    }
  }
  // Loop over locations and levels
  for (auto i = 0; i < nLocations_; i++) {
    for (auto k = 0; k < nLevels_; k++) {
      // Compute vectors of coeffs at each point i, k and store in coeffsFSet
      computeVectorsAt(i, k, linearEnsemble, linearErrors, coeffsFSet);
    }
  }
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HtlmCalculator<MODEL>::computeVectorsAt(const atlas::idx_t i,
                                             const atlas::idx_t k,
                                             const IncrementEnsemble_ & linearEnsemble,
                                             const IncrementEnsemble_ & linearErrors,
                                             atlas::FieldSet & coeffsFSet) const {
  // Make influenceMatrix at i, k
  EigenMatrix influenceMatrix = makeInfluenceMatrix(i, k, linearEnsemble);
  // Compute its singular value decomposition and obtain the U matrix of this
  const auto svd = Eigen::BDCSVD<EigenMatrix>(influenceMatrix * influenceMatrix.transpose(),
                                              Eigen::ComputeFullU | Eigen::ComputeFullV);
  const auto U = svd.matrixU();
  for (const auto & var : updateVars_.variables()) {
    // Copy (RMS-by-level scaled?) linearError for var at i, k into an EigenVector
    EigenVector linearErrorVector(ensembleSize_);
    for (auto m = 0; m < ensembleSize_; m++) {
      if (params_.rmsScaling) {  // apply RMS-by-level scaling
        linearErrorVector(m, 0)
          = atlas::array::make_view<double, 2>(linearErrors[m].fieldSet()[var])(i, k)
            / rmsVals_.at(var)[k];
      } else {  // don't apply RMS-by-level scaling
        linearErrorVector(m, 0)
          = atlas::array::make_view<double, 2>(linearErrors[m].fieldSet()[var])(i, k);
      }
    }
    // Compute vector of coeffs for var at i, k
    const EigenVector coeffs = computeVector(svd, U, influenceMatrix, linearErrorVector,
                                             regularization_->getRegularizationValue(var, i, k));
    // Copy coeffs into FieldSet
    auto coeffsFieldView = atlas::array::make_view<double, 3>(coeffsFSet[var]);
    for (auto x = 0; x < vectorSize_; x++) {
      coeffsFieldView(i, k, x) = coeffs[x];
    }
  }
}

//------------------------------------------------------------------------------

template<typename MODEL>
const typename HtlmCalculator<MODEL>::EigenMatrix HtlmCalculator<MODEL>::makeInfluenceMatrix(
                                                  const atlas::idx_t i,
                                                  const atlas::idx_t k,
                                                  const IncrementEnsemble_ & linearEnsemble) const {
  EigenMatrix influenceMatrix(vectorSize_, ensembleSize_);
  for (auto m = 0; m < ensembleSize_; m++) {
    for (size_t v = 0; v < updateVars_.size(); v++) {
      const auto linearEnsembleArray
        = atlas::array::make_view<double, 2>(linearEnsemble[m].fieldSet()[updateVars_[v]]);
      if (k >= halfInfluenceSize_ && k < nLevels_ - halfInfluenceSize_) {  // general case
        for (auto s = 0; s < influenceSize_; s++) {
          influenceMatrix(v * influenceSize_ + s, m)
            = linearEnsembleArray(i, k - halfInfluenceSize_ + s);
        }
      } else if (k < halfInfluenceSize_) {  // bottom of model
        for (auto s = 0; s < influenceSize_; s++) {
          influenceMatrix(v * influenceSize_ + s, m) = linearEnsembleArray(i, s);
        }
      } else {  // top of model
        for (auto s = 0; s < influenceSize_; s++) {
          influenceMatrix(v * influenceSize_ + s, m)
            = linearEnsembleArray(i, nLevels_ - influenceSize_ + s);
        }
      }
    }
  }
  return influenceMatrix;
}

//------------------------------------------------------------------------------

template<typename MODEL>
const typename HtlmCalculator<MODEL>::EigenVector HtlmCalculator<MODEL>::computeVector(
                                                         const Eigen::BDCSVD<EigenMatrix> & svd,
                                                         const EigenMatrix & U,
                                                         const EigenMatrix & influenceMatrix,
                                                         const EigenVector & linearErrorVector,
                                                         const double & regularizationValue) const {
  // Equation 26 in https://doi.org/10.1175/MWR-D-20-0088.1
  const EigenMatrix sigma
    = (svd.singularValues() + EigenVector::Constant(vectorSize_, regularizationValue)).asDiagonal();
  return U * sigma.completeOrthogonalDecomposition().pseudoInverse() * U.transpose()
           * influenceMatrix * linearErrorVector;
}

}  // namespace oops

#endif  // OOPS_GENERIC_HTLMCALCULATOR_H_
