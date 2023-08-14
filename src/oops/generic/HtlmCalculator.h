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
#include <vector>

#include "oops/base/IncrementEnsemble.h"
#include "oops/generic/HtlmRegularization.h"

namespace oops {

template <typename MODEL>
class HtlmCalculatorParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(HtlmCalculatorParameters, Parameters);

 public:
  /* The influence region size is the number of points that go into the calculation of
  coefficient vector for a tlm grid point. The influence region itself includes the tlm
  grid point and points above and below it in the vertical, evenly split when possible.
  For a grid point at say the bottom level, all influence points would be above it.
  The influence region also includes other model variables at the same locations.  */ 
  Parameter<HtlmRegularizationParameters> regularizationParams{"regularization", {}, this};
  Parameter<bool>                 rms{"rms scaling", true, this};
};

template <typename MODEL>
class HtlmCalculator{
  typedef HtlmCalculatorParameters<MODEL>   HtlmCalculatorParameters_;
  typedef Increment<MODEL>                  Increment_;
  typedef IncrementEnsemble<MODEL>          IncrementEnsemble_;
  typedef Geometry<MODEL>                   Geometry_;

 public:
  static const std::string classname() {return "oops::HtlmCalculator";}
  HtlmCalculator(const HtlmCalculatorParameters_ &, const Variables &, const size_t,
                 const Geometry_ &, const size_t, const size_t, const util::DateTime &,
                 const atlas::idx_t);
  void calcCoeffs(const IncrementEnsemble_ &,
                  const IncrementEnsemble_ &,
                  atlas::FieldSet &,
                  const atlas::FunctionSpace &) const;

 private:
  const HtlmCalculatorParameters_ params_;
  atlas::idx_t ensembleSize_;
  atlas::idx_t influenceSize_;
  atlas::idx_t halfInfluenceSize_;
  Variables vars_;
  std::unique_ptr<HtlmRegularization> regularizationPtr_;
  atlas::idx_t horizExt_;
  atlas::idx_t vertExt_;
};

//----------------------------------------------------------------------------------

template <typename MODEL>
HtlmCalculator<MODEL>::HtlmCalculator(const HtlmCalculatorParameters_ & params,
                                      const Variables & vars,
                                      const size_t ensembleSize,
                                      const Geometry_ & geomTLM,
                                      const size_t nLocations,
                                      const size_t nLevels,
                                      const util::DateTime & startTime,
                                      const atlas::idx_t influenceRegionSize)
: params_(params), ensembleSize_(ensembleSize), influenceSize_(influenceRegionSize),
  halfInfluenceSize_(influenceSize_/2), vars_(vars), horizExt_(nLocations), vertExt_(nLevels) {
  if (params_.regularizationParams.value().parts.value() == boost::none) {
    regularizationPtr_ = std::make_unique<HtlmRegularization>(params_.regularizationParams.value());
  } else {
    Increment_ regularizationIncrement(geomTLM, vars_, startTime);
    atlas::FieldSet regularizationFieldSet = regularizationIncrement.fieldSet();
    regularizationPtr_ =
      std::make_unique<HtlmRegularizationComponentDependent>(params_.regularizationParams.value(),
                                                             regularizationFieldSet);
  }
}

//-----------------------------------------------------------------------------------

// Runs the coefficient calculation looping over variables and grid points,
// also sizes the field set to store coefficient vectors and copies them.
template<typename MODEL>
void HtlmCalculator<MODEL>::calcCoeffs(const IncrementEnsemble_ & linearEnsemble,
                                       const IncrementEnsemble_ & linearErrorDe,
                                       atlas::FieldSet & coeffFieldSet,
                                       const atlas::FunctionSpace & fSpace) const {
  Log::trace() << "HtlmCalculator<MODEL>::coeffCalc() starting" << std::endl;
  // For each variable loop over every grid point and calculate the coefficient vector for each
  for (size_t varInd = 0; varInd < vars_.size(); ++varInd) {
     // make field set with size to store coefficient vectors
    for (size_t v = 0; v < vars_.size() * influenceSize_; v++) {
      atlas::Field coeffField = fSpace.createField<double>(
        atlas::option::name((vars_[varInd] + std::to_string(v)))
        | atlas::option::levels(vertExt_));
      coeffFieldSet.add(coeffField);
    }
    // get rms by level scaling
      const std::vector<double> rmsVals =
        params_.rms ? linearEnsemble[0].rmsByLevel(vars_[varInd]) : std::vector<double>{};
    // calculate coefficient vector for each grid point
    for (atlas::idx_t i = 0; i < horizExt_; ++i) {
      for (atlas::idx_t k = 0; k < vertExt_; ++k) {
        // matrices used for calculation of coeffs and storage
        // a set for each grid point for future parallelization
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>
                 influenceMat(influenceSize_*vars_.size(), ensembleSize_);
        Eigen::Matrix<double, Eigen::Dynamic, 1> linErrVec(ensembleSize_);
        Eigen::Matrix<double, Eigen::Dynamic, 1> coeffVect(influenceSize_*vars_.size());

        // Populate influenceMat (M) and linErrVec (delta e) declared in class
        Log::trace() << "HtlmCalculator<MODEL>::coeffCalc() filling matrices starting"
                                                                            << std::endl;
        for (atlas::idx_t ensInd = 0; ensInd < ensembleSize_; ++ensInd) {
          // populate linErrVec
          const atlas::FieldSet & linearEnsembleFset = linearEnsemble[ensInd].fieldSet();
          const atlas::FieldSet & linErrFset = linearErrorDe[ensInd].fieldSet();
          const auto linErrView = atlas::array::make_view<double, 2>(linErrFset[vars_[varInd]]);
          if (params_.rms) {
            linErrVec(ensInd, 0) = linErrView(i, k)/rmsVals[k];
          } else {
            linErrVec(ensInd, 0) = linErrView(i, k);
          }
          // populate influenceMat
          for (size_t varInd2 = 0; varInd2 < vars_.size(); ++varInd2) {
            const auto linearEnsembleView =
              atlas::array::make_view<double, 2>(linearEnsembleFset[vars_[varInd2]]);
            for (atlas::idx_t infInd = 0; infInd < influenceSize_; ++infInd) {
              if (k-halfInfluenceSize_ > 0 &&
                      k + halfInfluenceSize_ < linearEnsembleView.shape(1) ) {
                // middle case
                influenceMat(influenceSize_*varInd2 + infInd, ensInd) =
                        linearEnsembleView(i, k-halfInfluenceSize_ + infInd);
              } else if (k-halfInfluenceSize_ <= 0) {
                // start of increment edge case
                influenceMat(influenceSize_*varInd2 + infInd, ensInd) =
                                                    linearEnsembleView(i, infInd);
              } else if (k+halfInfluenceSize_ >= linearEnsembleView.shape(1)) {
                // end of increment edge case
                influenceMat(influenceSize_*varInd2 + infInd, ensInd) =
                    linearEnsembleView(i, (linearEnsembleView.shape(1)-influenceSize_) + infInd);
              } else {
                // Checking cases are correct
                Log::info() << "HtlmCalculator<MODEL>::coeffCalc() "
                              "unable to get influence region" << std::endl;
                abort();
              }  // end else
            }   // end infInd
          }   // end varInd2
        }   // end ensInd
        Log::trace() << "HtlmCalculator<MODEL>::coeffCalc() filling matrices done" << std::endl;

        // Calculate the coefficient vector for grid point at i,k
        Log::trace() << "HtlmCalculator<MODEL>::coeffCalc() calculating coefficients starting"
                                                                                    << std::endl;
        // Calculate the coefficient vector coeffVect for the grid point at i,k using
        // influenceMat (M) and linErrVec (delta e) declared in class
        const auto normalMat = influenceMat*influenceMat.transpose();
        const auto bdc_svd_solver =
                   Eigen::BDCSVD<Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic>>
                                          (normalMat, Eigen::ComputeFullU | Eigen::ComputeFullV);
        Eigen::Matrix<double, Eigen::Dynamic, Eigen::Dynamic> singularValues =
                                                  bdc_svd_solver.singularValues().asDiagonal();
        const auto U = bdc_svd_solver.matrixU();
        // Apply regularization
        for (int ij = 0; ij < singularValues.cols(); ++ij) {
          singularValues(ij, ij) += regularizationPtr_->getRegularizationValue(vars_[varInd], i, k);
        }
        const auto pseudoInv = singularValues.completeOrthogonalDecomposition().pseudoInverse();
        coeffVect = U*pseudoInv*U.transpose()*influenceMat*linErrVec;
        Log::trace() << "HtlmCalculator<MODEL>::coeffCalc() calculating coefficients done"
                                                                              << std::endl;

        // Copy the coeff vect into the its field set.
        Log::trace() << "HtlmCalculator<MODEL>::coeffCalc() placing vector starting" << std::endl;
        for (atlas::idx_t coeffInd = 0; coeffInd < atlas::idx_t(coeffVect.size()); ++coeffInd) {
          auto coeffsView = atlas::array::make_view<double, 2>(
            coeffFieldSet[vars_[varInd] + std::to_string(coeffInd)]);
          coeffsView(i, k) = coeffVect[coeffInd];
        }

        Log::trace() << "HtlmCalculator<MODEL>::coeffCalc() placing vector done" << std::endl;
      }  //  end for k
    }  //  end for i
  }  //  end for varInd
  Log::trace() << "HtlmCalculator<MODEL>::coeffCalc() done" << std::endl;
}

//-----------------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_GENERIC_HTLMCALCULATOR_H_
