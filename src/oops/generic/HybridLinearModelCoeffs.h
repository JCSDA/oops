/*
 * (C) Copyright 2022-2023 UCAR.
 * (C) Crown copyright 2022-2023 Met Office.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_GENERIC_HYBRIDLINEARMODELCOEFFS_H_
#define OOPS_GENERIC_HYBRIDLINEARMODELCOEFFS_H_

#include <map>
#include <vector>

#include "oops/generic/HtlmCalculator.h"
#include "oops/generic/HtlmEnsemble.h"

namespace oops {

template <typename MODEL>
class HybridLinearModelCoeffsParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(HybridLinearModelCoeffsParameters, Parameters);
  typedef HtlmCalculatorParameters<MODEL>    CalculatorParameters_;
  typedef HtlmEnsembleParameters<MODEL>      EnsembleParameters_;

 public:
  RequiredParameter<Variables> trainingVars{"training variables", this};
  RequiredParameter<atlas::idx_t> influenceSize{"influence region size", this};
  RequiredParameter<util::DateTime> windowBegin{"window begin", this};
  RequiredParameter<util::Duration> windowLength{"window length", this};
  OptionalParameter<EnsembleParameters_> ensemble{"ensemble", this};
  OptionalParameter<CalculatorParameters_> calculator{"calculator", this};
  OptionalParameter<bool> fromFile{"from file", this};  // TO DO
};

//------------------------------------------------------------------------------

template <typename MODEL>
class HybridLinearModelCoeffs {
 public:
  typedef Geometry<MODEL>                             Geometry_;
  typedef HtlmCalculator<MODEL>                       HtlmCalculator_;
  typedef HtlmEnsemble<MODEL>                         HtlmEnsemble_;
  typedef HtlmSimplifiedLinearModel<MODEL>            HtlmSimplifiedLinearModel_;
  typedef HybridLinearModelCoeffsParameters<MODEL>    Parameters_;
  typedef Increment<MODEL>                            Increment_;

  HybridLinearModelCoeffs(const Parameters_ &,
                          const Geometry_ &,
                          const util::Duration &,
                          HtlmSimplifiedLinearModel_ &);
  void updateIncTL(Increment_ &) const;
  void updateIncAD(Increment_ &) const;

 private:
  void generate(const Geometry_ &, const util::Duration &, HtlmSimplifiedLinearModel_ &);

  const Parameters_ & params_;
  const Variables trainingVars_;
  const atlas::idx_t nLevels_;
  const atlas::idx_t influenceSize_;
  atlas::Field updateStencil_;
  std::map<util::DateTime, atlas::FieldSet> coeffsSaver_;
};

//------------------------------------------------------------------------------

template<typename MODEL>
HybridLinearModelCoeffs<MODEL>::HybridLinearModelCoeffs(
                                                 const Parameters_ & params,
                                                 const Geometry_ & updateGeometry,
                                                 const util::Duration & updateTstep,
                                                 HtlmSimplifiedLinearModel_ & simplifiedLinearModel)
: params_(params), trainingVars_(params_.trainingVars),
  nLevels_(updateGeometry.variableSizes(trainingVars_)[0]), influenceSize_(params_.influenceSize),
  updateStencil_("update stencil", atlas::array::make_datatype<int>(),
                 atlas::array::make_shape(nLevels_, influenceSize_)) {
  if (updateTstep % simplifiedLinearModel.timeResolution() != 0) {
    ABORT("HybridLinearModelCoeffs<MODEL>::HybridLinearModelCoeffs: "
          "update tstep is not a multiple of simplified linear model tstep");
  }
  // Make stencil for applying coefficients
  auto stencilView = atlas::array::make_view<atlas::idx_t, 2>(updateStencil_);
  const atlas::idx_t halfInfluenceSize = influenceSize_ / 2;
  for (atlas::idx_t k = 0; k < nLevels_; k++) {
    for (atlas::idx_t infInd = 0; infInd < influenceSize_; infInd++) {
      if (k - halfInfluenceSize > 0 && k + halfInfluenceSize < nLevels_) {  // Middle case
        stencilView(k, infInd) = k - halfInfluenceSize + infInd;
      } else if (k - halfInfluenceSize <= 0) {  // Bottom case
        stencilView(k, infInd) = infInd;
      } else if (k + halfInfluenceSize >= nLevels_) {  // Top case
        stencilView(k, infInd) = (nLevels_ - influenceSize_) + infInd;
      } else {
        ABORT("HybridLinearModelCoeffs<MODEL>::HybridLinearModelCoeffs(): "
              "unable to determine position in influence region");
      }
    }
  }
  // Determine source of and obtain coefficients
  if (params_.ensemble.value() != boost::none && params_.calculator.value() != boost::none) {
    generate(updateGeometry, updateTstep, simplifiedLinearModel);
  } else if (params.fromFile.value() != boost::none) {
    // TO DO
  } else {
    ABORT("HybridLinearModelCoeffs<MODEL>::HybridLinearModelCoeffs(): no source of coefficients");
  }
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::generate(const Geometry_ & updateGeometry,
                                              const util::Duration & updateTstep,
                                              HtlmSimplifiedLinearModel_ & simplifiedLinearModel) {
  HtlmEnsemble_ ensemble(*params_.ensemble.value(), simplifiedLinearModel.variables(),
                         updateGeometry);
  HtlmCalculator_ calculator(*params_.calculator.value(), trainingVars_, ensemble.size(),
                             updateGeometry, updateGeometry.functionSpace().size(), nLevels_,
                             params_.windowBegin, influenceSize_);
  util::DateTime time(params_.windowBegin);
  while (time < (params_.windowBegin.value() + params_.windowLength.value())) {
    time += updateTstep;
    ensemble.step(updateTstep, simplifiedLinearModel);
    // Create an empty atlas::FieldSet to store coefficients at this time
    atlas::FieldSet coeffsFieldSet;
    // Store this in the map
    coeffsSaver_.emplace(time, coeffsFieldSet);
    // Calculate the coefficients and assign them to the FieldSet
    calculator.calcCoeffs(ensemble.getLinearEns(), ensemble.getLinearErrDe(), coeffsFieldSet);
    // Update increments with coefficients before next step
    for (size_t m = 0; m < ensemble.size(); m++) {
      updateIncTL(ensemble.getLinearEns()[m]);
    }
  }
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::updateIncTL(Increment_ & dx) const {
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::updateIncTL() starting" << std::endl;
  auto stencilView = atlas::array::make_view<int, 2>(updateStencil_);
  atlas::FieldSet & dxFset = dx.fieldSet();
  auto dxView = atlas::array::make_view<double, 2>(dxFset[trainingVars_[0]]);
  for (atlas::idx_t i = 0; i < dxView.shape(0); i++) {
    std::vector<double> updateVals(dxView.shape(1) * trainingVars_.size(), 0.0);
    // Calculate update values
    for (size_t varInd = 0; varInd < trainingVars_.size(); varInd++) {
      auto coeffView =
        atlas::array::make_view<double, 3>(coeffsSaver_.at(dx.validTime())[trainingVars_[varInd]]);
      for (atlas::idx_t k = 0; k < dxView.shape(1); k++) {
        for (size_t varInd2 = 0; varInd2 < trainingVars_.size(); varInd2++) {
          auto dxView = atlas::array::make_view<double, 2>(dxFset[trainingVars_[varInd2]]);
          for (atlas::idx_t coeffInd = 0; coeffInd < influenceSize_; coeffInd++) {
            updateVals[k + varInd * dxView.shape(1)] +=
              coeffView(i, k, varInd2 * influenceSize_ + coeffInd)
              * dxView(i, stencilView(k, coeffInd));
          }
        }
      }
    }
    // Update column
    for (size_t varInd = 0; varInd < trainingVars_.size(); varInd++) {
      auto dxView = atlas::array::make_view<double, 2>(dxFset[trainingVars_[varInd]]);
      for (atlas::idx_t k = 0; k < dxView.shape(1); k++) {
        dxView(i, k) += updateVals[k + varInd * dxView.shape(1)];
      }
    }
  }
  dx.synchronizeFields();
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::updateIncTL() done" << std::endl;
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::updateIncAD(Increment_ & dx) const {
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::updateIncAD() starting" << std::endl;
  auto stencilView = atlas::array::make_view<int, 2>(updateStencil_);
  atlas::FieldSet & dxFset = dx.fieldSet();
  auto dxView = atlas::array::make_view<double, 2>(dxFset[trainingVars_[0]]);
  for (atlas::idx_t i = 0; i < dxView.shape(0); i++) {
    std::vector<double> updateVals(dxView.shape(1) * trainingVars_.size(), 0.0);
    // Adjoint of "Update column"
    for (size_t varInd = 0; varInd < trainingVars_.size(); varInd++) {
      auto dxView = atlas::array::make_view<double, 2>(dxFset[trainingVars_[varInd]]);
      for (atlas::idx_t k = 0; k < dxView.shape(1); k++) {
        updateVals[k + varInd * dxView.shape(1)] += dxView(i, k);
      }
    }
    // Adjoint of "Calculate update values"
    for (size_t varInd = 0; varInd < trainingVars_.size(); ++varInd) {
      auto coeffView =
        atlas::array::make_view<double, 3>(coeffsSaver_.at(dx.validTime())[trainingVars_[varInd]]);
      for (atlas::idx_t k = 0; k < dxView.shape(1); k++) {
        for (size_t varInd2 = 0; varInd2 < trainingVars_.size(); varInd2++) {
          auto dxView = atlas::array::make_view<double, 2>(dxFset[trainingVars_[varInd2]]);
          for (atlas::idx_t coeffInd = 0; coeffInd < influenceSize_; coeffInd++) {
            dxView(i, stencilView(k, coeffInd)) +=
              coeffView(i, k, varInd2 * influenceSize_ + coeffInd)
              * updateVals[k + varInd * dxView.shape(1)];
          }
        }
      }
    }
  }
  dx.synchronizeFields();
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::updateIncAD() done" << std::endl;
}

}  // namespace oops

#endif  // OOPS_GENERIC_HYBRIDLINEARMODELCOEFFS_H_
