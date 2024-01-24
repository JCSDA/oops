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
#include <memory>
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
  RequiredParameter<Variables> updateVars{"update variables", this};
  RequiredParameter<atlas::idx_t> influenceSize{"influence region size", this};
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
  typedef HybridLinearModelCoeffsParameters<MODEL>    Parameters_;
  typedef Increment<MODEL>                            Increment_;
  typedef SimpleLinearModel<MODEL>                    SimpleLinearModel_;

  HybridLinearModelCoeffs(const eckit::Configuration &, const Geometry_ &, const util::Duration &);
  void obtain(SimpleLinearModel_ &, const Variables &);
  void updateIncTL(Increment_ &) const;
  void updateIncAD(Increment_ &) const;

 private:
  void makeCoeffsSaver();
  void makeUpdateStencil();
  void generate(SimpleLinearModel_ &, const Variables &);
  void read();
  void write() const;

  Parameters_ params_;
  const Variables updateVars_;
  const Geometry_ & updateGeometry_;
  const util::Duration & updateTstep_;
  const atlas::idx_t nLocations_;
  const atlas::idx_t nLevels_;
  const atlas::idx_t influenceSize_;
  atlas::Field updateStencil_;
  std::map<util::DateTime, atlas::FieldSet> coeffsSaver_;
  const util::TimeWindow timeWindow_;
};

//------------------------------------------------------------------------------

template<typename MODEL>
HybridLinearModelCoeffs<MODEL>::HybridLinearModelCoeffs(
                                                 const eckit::Configuration & config,
                                                 const Geometry_ & updateGeometry,
                                                 const util::Duration & updateTstep)
: params_(), updateVars_(config, "update variables"), updateGeometry_(updateGeometry),
  updateTstep_(updateTstep), nLocations_(updateGeometry_.functionSpace().size()),
  nLevels_(updateGeometry.variableSizes(updateVars_)[0]),
  influenceSize_(config.getInt("influence region size")),
  updateStencil_("update stencil", atlas::array::make_datatype<int>(),
                 atlas::array::make_shape(nLevels_, influenceSize_)),
  timeWindow_(config.getSubConfiguration("time window"))
{
  if (influenceSize_ % 2 == 0) {
    oops::Log::warning() << "HybridLinearModelCoeffs<MODEL>::HybridLinearModelCoeffs: "
                            "influence region size is not an odd number;"
                            "influence regions will not be centred on point of interest";
  }
  params_.deserialize(config);
  // Set up storage for coefficients
  makeCoeffsSaver();
  // Set up stencil for applying coefficients
  makeUpdateStencil();
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::obtain(SimpleLinearModel_ & simpleLinearModel,
                                            const Variables & vars) {
  // Determine source of and obtain coefficients
  if (params_.ensemble.value() != boost::none && params_.calculator.value() != boost::none) {
    generate(simpleLinearModel, vars);
  } else if (params_.input.value() != boost::none) {
    read();
  } else {
    ABORT("HybridLinearModelCoeffs<MODEL>::obtain(): no source of coefficients");
  }
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::makeCoeffsSaver() {
  // Create Fields for coeffs at each time, using FunctionSpace from updateGeometry
  util::DateTime time(timeWindow_.start());
  while (time < timeWindow_.end()) {
    time += updateTstep_;
    atlas::FieldSet coeffsFSet;
    coeffsSaver_.emplace(time, coeffsFSet);
    for (const auto & var : updateVars_.variables()) {
      coeffsSaver_.at(time).add(updateGeometry_.functionSpace().template createField<double>(
        atlas::option::name(var) | atlas::option::levels(nLevels_) |
        atlas::option::vector(influenceSize_ * updateVars_.size())));
    }
  }
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::makeUpdateStencil() {
  auto updateStencilArray = atlas::array::make_view<atlas::idx_t, 2>(updateStencil_);
  const auto halfInfluenceSize = influenceSize_ / 2;
  for (auto s = 0; s < influenceSize_; s++) {
    for (auto k = 0; k < halfInfluenceSize; k++) {  // bottom of model
      updateStencilArray(k, s) = s;
    }
    for (auto k = halfInfluenceSize; k < nLevels_ - halfInfluenceSize; k++) {  // general case
      updateStencilArray(k, s) = k - halfInfluenceSize + s;
    }
    for (auto k = nLevels_ - halfInfluenceSize; k < nLevels_; k++) {  // top of model
      updateStencilArray(k, s) = nLevels_ - influenceSize_ + s;
    }
  }
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::generate(SimpleLinearModel_ & simpleLinearModel,
                                              const Variables & vars) {
  HtlmEnsemble_ ensemble(*params_.ensemble.value(), simpleLinearModel, updateGeometry_, vars);
  HtlmCalculator_ calculator(*params_.calculator.value(), updateVars_, updateGeometry_,
                             influenceSize_, ensemble.size());
  util::DateTime time(timeWindow_.start());
  while (time < timeWindow_.end()) {
    time += updateTstep_;
    ensemble.step(updateTstep_, simpleLinearModel);
    calculator.setOfCoeffs(ensemble.getLinearEnsemble(), ensemble.getLinearErrors(),
                           coeffsSaver_.at(time));
    // Update increments with coefficients before next step
    for (size_t m = 0; m < ensemble.size(); m++) {
      updateIncTL(ensemble.getLinearEnsemble()[m]);
    }
  }
  if (params_.output.value() != boost::none) {
    write();
  }
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::read() {
  const std::vector<size_t> nLevelsAll(updateVars_.size(), nLevels_);
  eckit::LocalConfiguration inputConfig(*params_.input.value());
  const std::string baseFilepath = inputConfig.getString("base filepath");
  util::DateTime time(timeWindow_.start());
  while (time < timeWindow_.end()) {
    time += updateTstep_;
    const std::string filepath = baseFilepath + "_" + time.toStringIO();
    inputConfig.set("filepath", filepath);
    util::readFieldSet(updateGeometry_.getComm(), updateGeometry_.functionSpace(), nLevelsAll,
                       updateVars_.variables(), inputConfig, coeffsSaver_.at(time));
  }
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::write() const {
  eckit::LocalConfiguration outputConfig(*params_.output.value());
  const std::string baseFilepath = outputConfig.getString("base filepath");
  for (const auto & element : coeffsSaver_) {
    const std::string filepath = baseFilepath + "_" + element.first.toStringIO();
    outputConfig.set("filepath", filepath);
    util::writeFieldSet(updateGeometry_.getComm(), outputConfig, element.second);
  }
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::updateIncTL(Increment_ & dx) const {
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::updateIncTL() starting" << std::endl;
  const auto updateStencilArray = atlas::array::make_view<int, 2>(updateStencil_);
  atlas::FieldSet & dxFSet = dx.fieldSet().fieldSet();
  std::vector<double> updateVals(nLevels_ * updateVars_.size());
  for (auto i = 0; i < nLocations_; i++) {
    std::fill(updateVals.begin(), updateVals.end(), 0.0);
    // Calculate update values
    for (size_t v = 0; v < updateVars_.size(); v++) {
      auto coeffsView
        = atlas::array::make_view<double, 3>(coeffsSaver_.at(dx.validTime())[updateVars_[v]]);
      for (auto k = 0; k < nLevels_; k++) {
        for (size_t v2 = 0; v2 < updateVars_.size(); v2++) {
          auto dxArray = atlas::array::make_view<double, 2>(dxFSet[updateVars_[v2]]);
          for (auto s = 0; s < influenceSize_; s++) {
            updateVals[k + v * nLevels_]
              += coeffsView(i, k, v2 * influenceSize_ + s) * dxArray(i, updateStencilArray(k, s));
          }
        }
      }
    }
    // Update column
    for (size_t v = 0; v < updateVars_.size(); v++) {
      auto dxArray = atlas::array::make_view<double, 2>(dxFSet[updateVars_[v]]);
      for (auto k = 0; k < nLevels_; k++) {
        dxArray(i, k) += updateVals[k + v * nLevels_];
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
  const auto updateStencilArray = atlas::array::make_view<int, 2>(updateStencil_);
  atlas::FieldSet & dxFSet = dx.fieldSet().fieldSet();
  std::vector<double> updateVals(nLevels_ * updateVars_.size());
  for (auto i = 0; i < nLocations_; i++) {
    std::fill(updateVals.begin(), updateVals.end(), 0.0);
    // Adjoint of "Update column"
    for (size_t v = 0; v < updateVars_.size(); v++) {
      auto dxArray = atlas::array::make_view<double, 2>(dxFSet[updateVars_[v]]);
      for (auto k = 0; k < nLevels_; k++) {
        updateVals[k + v * nLevels_] += dxArray(i, k);
      }
    }
    // Adjoint of "Calculate update values"
    for (size_t v = 0; v < updateVars_.size(); v++) {
      auto coeffsView
        = atlas::array::make_view<double, 3>(coeffsSaver_.at(dx.validTime())[updateVars_[v]]);
      for (auto k = 0; k < nLevels_; k++) {
        for (size_t v2 = 0; v2 < updateVars_.size(); v2++) {
          auto dxArray = atlas::array::make_view<double, 2>(dxFSet[updateVars_[v2]]);
          for (auto s = 0; s < influenceSize_; s++) {
            dxArray(i, updateStencilArray(k, s))
              += coeffsView(i, k, v2 * influenceSize_ + s) * updateVals[k + v * nLevels_];
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
