/*
 * (C) Copyright 2022 MetOffice.
 * (C) Copyright 2021-2022 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#ifndef OOPS_GENERIC_HTLMCALCULATOR_H_
#define OOPS_GENERIC_HTLMCALCULATOR_H_

#include <Eigen/Geometry>
#include <Eigen/SVD>
#include <memory>
#include <string>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/Variables.h"
#include "oops/generic/HtlmRegularization.h"
#include "oops/mpi/mpi.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace oops {

// Parameters for coefficient calculation
template <typename MODEL>
class HtlmCalculatorParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(HtlmCalculatorParameters, Parameters);
  typedef typename Geometry<MODEL>::Parameters_ GeometryParameters_;

 public:
  /* The influence region size is the number of points that go into the calculation of
  coefficient vector for a tlm grid point. The influence region itself includes the tlm
  grid point and points above and below it in the vertical, evenly split when possible.
  For a grid point at say the bottom level, all influence points would be above it.
  The influence region also includes other model variables at the same locations.  */ 
  RequiredParameter<atlas::idx_t> influenceRegionSize{"influence region size", this};
  Parameter<HtlmRegularizationParameters> regularizationParams{"regularization", {}, this};
  Parameter<bool>                 rms{"rms scaling", true, this};
};

// Class declarations for HtlmCalculator
template <typename MODEL>
class HtlmCalculator{
  typedef HtlmCalculatorParameters<MODEL>   HtlmCalculatorParameters_;
  typedef Increment<MODEL>                  Increment_;
  typedef IncrementEnsemble<MODEL>          IncrementEnsemble_;
  typedef Geometry<MODEL>                   Geometry_;

 public:
  static const std::string classname() {return "oops::HtlmCalculator";}

  HtlmCalculator(const HtlmCalculatorParameters_ &, const Variables &, const size_t &,
                                                     const Geometry_ &, const util::DateTime &);

  void calcCoeffs(const IncrementEnsemble_ &,
                 const IncrementEnsemble_ &,
                 atlas::FieldSet &);

 private:
  const HtlmCalculatorParameters_ params_;
  atlas::idx_t ensembleSize_;
  atlas::idx_t influenceSize_;
  atlas::idx_t halfInfluenceSize_;
  Variables vars_;

  // regularization paramater in coeff calculation
  std::unique_ptr<HtlmRegularization> regularizationPtr_;

  // Increment geometry used to get needed dimensions.
  const Geometry_ & incrementGeometry_;
  atlas::idx_t horizExt_;
  atlas::idx_t vertExt_;
};

//----------------------------------------------------------------------------------

template <typename MODEL>
HtlmCalculator<MODEL>::HtlmCalculator(const HtlmCalculatorParameters_ & params,
                                                      const Variables & vars,
                                                      const size_t & ensembleSize,
                                                      const Geometry_ & geomTLM,
                                                      const util::DateTime & startTime)
: params_(params), ensembleSize_(atlas::idx_t(ensembleSize)),
influenceSize_(params_.influenceRegionSize.value()), halfInfluenceSize_(influenceSize_/2),
vars_(vars),
 incrementGeometry_(geomTLM) {
  // Use increment geometry to get needed dimensions
  Increment_ coeffSetupInc(incrementGeometry_, vars_, startTime);
  atlas::FieldSet coeffSetupFset = coeffSetupInc.fieldSet();
  horizExt_ = coeffSetupFset[0].shape(0);
  vertExt_ = coeffSetupFset[0].shape(1);

  // Set up regularization
  if (params_.regularizationParams.value().parts.value() == boost::none) {
    regularizationPtr_ = std::make_unique<HtlmRegularization>(params_.regularizationParams.value());
  } else {
    regularizationPtr_ =
      std::make_unique<HtlmRegularizationComponentDependent>(params_.regularizationParams.value(),
                                                             coeffSetupFset);
  }
}

//-----------------------------------------------------------------------------------

// Runs the coefficient calculation looping over variables and grid points,
// also sizes the field set to store coefficient vectors and copies them.
template<typename MODEL>
void HtlmCalculator<MODEL>::calcCoeffs(const IncrementEnsemble_ & linearEnsemble,
                                       const IncrementEnsemble_ & linearErrorDe,
                                       atlas::FieldSet & coeffFieldSet) {
  Log::trace() << "HtlmCalculator<MODEL>::coeffCalc() starting" << std::endl;
  // For each variable loop over every grid point and calculate the coefficient vector for each
  for (size_t varInd = 0; varInd < vars_.size(); ++varInd) {
     // make field set with size to store coefficient vectors
      atlas::Field
         coeffField(vars_[varInd], atlas::array::make_datatype<double>(),
            atlas::array::make_shape(horizExt_, vertExt_, vars_.size() * influenceSize_));
      coeffFieldSet.add(coeffField);
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
                influenceMat(vars_.size()*varInd2 + infInd, ensInd) =
                        linearEnsembleView(i, k-halfInfluenceSize_ + infInd);
              } else if (k-halfInfluenceSize_ <= 0) {
                // start of increment edge case
                influenceMat(vars_.size()*varInd2 + infInd, ensInd) =
                                                    linearEnsembleView(i, infInd);
              } else if (k+halfInfluenceSize_ >= linearEnsembleView.shape(1)) {
                // end of increment edge case
                influenceMat(vars_.size()*varInd2 + infInd, ensInd) =
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
        auto coeffsView = atlas::array::make_view<double, 3>(coeffFieldSet[vars_[varInd]]);
        for (atlas::idx_t coeffInd = 0; coeffInd < atlas::idx_t(coeffVect.size()); ++coeffInd) {
          coeffsView(i, k, coeffInd) = coeffVect[coeffInd];
        }
        Log::trace() << "HtlmCalculator<MODEL>::coeffCalc() placing vector done" << std::endl;
      }  //  end for k
    }  //  end for i
  }  //  end for varInd
  Log::trace() << "HtlmCalculator<MODEL>::coeffCalc() done" << std::endl;
}

//-----------------------------------------------------------------------------------

}  //  namespace oops
#endif  //  OOPS_GENERIC_HTLMCALCULATOR_H_
