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
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/ParameterTraitsVariables.h"
#include "oops/base/Variables.h"
#include "oops/generic/HtlmCalculator.h"
#include "oops/generic/HtlmEnsemble.h"
#include "oops/util/FieldSetHelpers.h"

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
  OptionalParameter<eckit::LocalConfiguration> input{"input", this};
  OptionalParameter<eckit::LocalConfiguration> output{"output", this};
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

  HybridLinearModelCoeffs(const eckit::Configuration &, const Geometry_ &,
                          const util::Duration &, HtlmSimplifiedLinearModel_ &);
  void updateIncTL(Increment_ &) const;
  void updateIncAD(Increment_ &) const;

 private:
  void generate(const Geometry_ &, const util::Duration &, HtlmSimplifiedLinearModel_ &);
  void read(const Geometry_ &, const util::Duration &);
  void write(const eckit::mpi::Comm &) const;

  Parameters_ params_;
  const Variables trainingVars_;
  const atlas::idx_t nLevels_;
  const atlas::idx_t influenceSize_;
  atlas::Field updateStencil_;
  std::map<util::DateTime, atlas::FieldSet> coeffsSaver_;
};

//------------------------------------------------------------------------------

template<typename MODEL>
HybridLinearModelCoeffs<MODEL>::HybridLinearModelCoeffs(const eckit::Configuration & config,
                                                        const Geometry_ & updateGeometry,
                                                        const util::Duration & updateTstep,
                                                 HtlmSimplifiedLinearModel_ & simplifiedLinearModel)
: params_(), trainingVars_(config, "training variables"),
  nLevels_(updateGeometry.variableSizes(trainingVars_)[0]),
  influenceSize_(config.getInt("influence region size")),
  updateStencil_("update stencil", atlas::array::make_datatype<int>(),
                 atlas::array::make_shape(nLevels_, influenceSize_))
{
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
  params_.deserialize(config);
  if (params_.ensemble.value() != boost::none && params_.calculator.value() != boost::none) {
    generate(updateGeometry, updateTstep, simplifiedLinearModel);
  } else if (params_.input.value() != boost::none) {
    read(updateGeometry, updateTstep);
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
    calculator.calcCoeffs(ensemble.getLinearEns(), ensemble.getLinearErrDe(), coeffsFieldSet,
                          updateGeometry.functionSpace());
    // Update increments with coefficients before next step
    for (size_t m = 0; m < ensemble.size(); m++) {
      updateIncTL(ensemble.getLinearEns()[m]);
    }
  }
  if (params_.output.value() != boost::none) {
    write(updateGeometry.getComm());
  }
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::read(const Geometry_ & updateGeometry,
                                          const util::Duration & updateTstep) {
  eckit::LocalConfiguration inputConfig(*params_.input.value());
  const std::string baseFilepath = inputConfig.getString("base filepath");
  Variables fileVars;
  for (const std::string & var : trainingVars_.variables()) {
    for (size_t v = 0; v < trainingVars_.size() * influenceSize_; v++) {
      fileVars.push_back(var + std::to_string(v));
      fileVars.addMetaData(var + std::to_string(v), "levels", nLevels_);
    }
  }
  util::DateTime time(params_.windowBegin);
  while (time < (params_.windowBegin.value() + params_.windowLength.value())) {
    time += updateTstep;
    const std::string filepath = baseFilepath + "_" + time.toStringIO();
    inputConfig.set("filepath", filepath);
    atlas::FieldSet coeffsFieldSet;
    coeffsSaver_.emplace(time, coeffsFieldSet);
    util::readFieldSet(updateGeometry.getComm(), updateGeometry.functionSpace(), fileVars,
                       inputConfig, coeffsFieldSet);
  }
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::write(const eckit::mpi::Comm & comm) const {
  eckit::LocalConfiguration outputConfig(*params_.output.value());
  const std::string baseFilepath = outputConfig.getString("base filepath");
  for (const auto & element : coeffsSaver_) {
    const std::string filepath = baseFilepath + "_" + element.first.toStringIO();
    outputConfig.set("filepath", filepath);
    util::writeFieldSet(comm, outputConfig, element.second);
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
      for (atlas::idx_t k = 0; k < dxView.shape(1); k++) {
        for (size_t varInd2 = 0; varInd2 < trainingVars_.size(); varInd2++) {
          auto dxView = atlas::array::make_view<double, 2>(dxFset[trainingVars_[varInd2]]);
          for (atlas::idx_t coeffInd = 0; coeffInd < influenceSize_; coeffInd++) {
            auto coeffView = atlas::array::make_view<double, 2>(coeffsSaver_.at(dx.validTime())[
              trainingVars_[varInd] + std::to_string(varInd2 * influenceSize_ + coeffInd)]);
            updateVals[k + varInd * dxView.shape(1)] +=
              coeffView(i, k)
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
      for (atlas::idx_t k = 0; k < dxView.shape(1); k++) {
        for (size_t varInd2 = 0; varInd2 < trainingVars_.size(); varInd2++) {
          auto dxView = atlas::array::make_view<double, 2>(dxFset[trainingVars_[varInd2]]);
          for (atlas::idx_t coeffInd = 0; coeffInd < influenceSize_; coeffInd++) {
            auto coeffView = atlas::array::make_view<double, 2>(coeffsSaver_.at(dx.validTime())[
              trainingVars_[varInd] + std::to_string(varInd2 * influenceSize_ + coeffInd)]);
            dxView(i, stencilView(k, coeffInd)) +=
              coeffView(i, k)
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
