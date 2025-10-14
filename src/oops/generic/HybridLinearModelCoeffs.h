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
#include "oops/base/Variables.h"
#include "oops/generic/HtlmCalculator.h"
#include "oops/util/FieldSetHelpers.h"
#include "oops/util/for_each.h"
#include "oops/util/ParallelFieldSetIO.h"
#include "oops/util/Timer.h"

namespace oops {

/*
 * Configuration options for HybridLinearModelCoeffs:
 *
 * This configuration defines the behavior of the hybrid linear model coefficients,
 * including the variables to update, the spatial influence region, and optional
 * components like ensemble settings and HTLM calculator parameters.
 *
 * Keys:
 * ─────────────────────────────────────────────────────────────────────────────
 * "update variables"     : (Required) List of variables to be updated by the model.
 *
 * "influence region size": (Required) Integer specifying the size of the spatial
 *                          region (in grid points or levels) that influences each update.
 *
 * "ensemble"             : (Optional) Configuration block for ensemble settings.
 *
 * "calculator"           : (Optional) Configuration block for HTLM calculator settings.
 *
 * "input"                : (Optional) Configuration block for input sources or overrides.
 *
 * "output"               : (Optional) Configuration block for output settings.
 */

//------------------------------------------------------------------------------

template <typename MODEL>
class HybridLinearModelCoeffs {
 public:
  typedef Geometry<MODEL>                             Geometry_;
  typedef HtlmCalculator<MODEL>                       HtlmCalculator_;
  typedef HtlmEnsemble<MODEL>                         HtlmEnsemble_;
  typedef Increment<MODEL>                            Increment_;
  typedef SimpleLinearModel<MODEL>                    SimpleLinearModel_;

  HybridLinearModelCoeffs(const eckit::Configuration &, const Geometry_ &, const util::Duration &);
  static const std::string classname() {return "oops::HybridLinearModelCoeffs";}
  void obtain(SimpleLinearModel_ &, const Variables &);
  void updateIncTL(Increment_ &) const;
  void updateIncAD(Increment_ &) const;

 private:
  void makeCoeffsSaver();
  void makeUpdateStencil();
  void generate(SimpleLinearModel_ &, const Variables &);
  void read();
  void write() const;

  const eckit::LocalConfiguration coeffsConfig_;
  const Variables updateVars_;
  const Geometry_ & updateGeometry_;
  const util::Duration & updateTstep_;
  const atlas::idx_t nLocations_;
  const atlas::idx_t nLevels_;
  const atlas::idx_t influenceSize_;
  atlas::Field updateStencil_;
  std::map<util::DateTime, atlas::FieldSet> coeffsSaver_;
  const util::TimeWindow timeWindow_;
  const std::vector<atlas::idx_t> owned_;
};

//------------------------------------------------------------------------------

template<typename MODEL>
HybridLinearModelCoeffs<MODEL>::HybridLinearModelCoeffs(
                                                 const eckit::Configuration & config,
                                                 const Geometry_ & updateGeometry,
                                                 const util::Duration & updateTstep)
: coeffsConfig_(config), updateVars_(config, "update variables"), updateGeometry_(updateGeometry),
  updateTstep_(updateTstep), nLocations_(updateGeometry_.functionSpace().size()),
  nLevels_(updateGeometry.variableSizes(updateVars_)[0]),
  influenceSize_(config.getInt("influence region size")),
  updateStencil_("update stencil", atlas::array::make_datatype<int>(),
                 atlas::array::make_shape(nLevels_, influenceSize_)),
  timeWindow_(config.getSubConfiguration("time window")),
  owned_([&]() {
    std::vector<atlas::idx_t> owned;
    const auto ownedView = atlas::array::make_view<int, 2>(updateGeometry_.fields()["owned"]);
    for (auto i = 0; i < nLocations_; i++) {
      if (ownedView(i, 0) > 0) owned.push_back(i);
    }
    return owned;
  }())
{
  if (influenceSize_ % 2 == 0) {
    ABORT("HybridLinearModelCoeffs<MODEL>::HybridLinearModelCoeffs: "
          "influence region size is not an odd number;"
          "influence regions will not be centred on point of interest");
  }
  // Set up storage for coefficients
  makeCoeffsSaver();
  // Set up stencil for applying coefficients
  makeUpdateStencil();
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::obtain(SimpleLinearModel_ & simpleLinearModel,
                                            const Variables & vars) {
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::obtain() starting" << std::endl;
  // Determine source of and obtain coefficients
  if (coeffsConfig_.has("ensemble") && coeffsConfig_.has("calculator")) {
    generate(simpleLinearModel, vars);
  } else if (coeffsConfig_.has("input")) {
    read();
  } else {
    ABORT("HybridLinearModelCoeffs<MODEL>::obtain(): no source of coefficients");
  }
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::obtain() done" << std::endl;
  if (coeffsConfig_.has("output")) {
    write();
  }
}


//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::makeCoeffsSaver() {
  // Create Fields for coeffs at each time, using FunctionSpace from updateGeometry
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::makeCoeffsSaver() starting" << std::endl;
  util::DateTime time(timeWindow_.start());
  while (time < timeWindow_.end()) {
    time += updateTstep_;
    atlas::FieldSet coeffsFSet;
    coeffsSaver_.emplace(time, coeffsFSet);
    for (const auto & var : updateVars_) {
      coeffsSaver_.at(time).add(updateGeometry_.functionSpace().template createField<double>(
        atlas::option::halo(0) | atlas::option::name(var.name()) | atlas::option::levels(nLevels_) |
        atlas::option::vector(influenceSize_ * updateVars_.size())));
    }
  }
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::makeCoeffsSaver() done" << std::endl;
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::makeUpdateStencil() {
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::makeUpdateStencil() starting" << std::endl;
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
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::makeUpdateStencil() done" << std::endl;
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::generate(SimpleLinearModel_ & simpleLinearModel,
                                              const Variables & vars) {
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::generate() starting" << std::endl;
  eckit::LocalConfiguration ensConf(coeffsConfig_, "ensemble");
  eckit::LocalConfiguration calcConf(coeffsConfig_, "calculator");
  HtlmEnsemble_ ensemble(ensConf,
                         simpleLinearModel, updateGeometry_, vars);
  HtlmCalculator_ calculator(calcConf, updateVars_,
                             updateGeometry_, influenceSize_, ensemble, owned_);
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
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::generate() done" << std::endl;
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::read() {
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::read() starting" << std::endl;
  util::Timer timer(classname(), "read");
  const std::vector<size_t> nLevelsAll(updateVars_.size(), nLevels_);
  eckit::LocalConfiguration inputConfig(coeffsConfig_, "input");
  const std::string baseFilepath = inputConfig.getString("base filepath");
  util::DateTime time(timeWindow_.start());
  if (inputConfig.getBool("legacy", true)) {
    while (time < timeWindow_.end()) {
      time += updateTstep_;
      const std::string filepath = baseFilepath + "_" + time.toStringIO();
      inputConfig.set("filepath", filepath);
      // TODO(someone): when updateVars_ have levels, replace by the call taking vars
      util::readFieldSet(updateGeometry_.getComm(), updateGeometry_.functionSpace(), nLevelsAll,
                         updateVars_.variables(), inputConfig, coeffsSaver_.at(time));
    }
  } else {
    util::ParallelFieldSetIO io(updateGeometry_.functionSpace(),
                                inputConfig.getString("grid name"),
                                util::ParallelFieldSetIO::Mode::Read);
    while (time < timeWindow_.end()) {
      time += updateTstep_;
      const std::string filepath = baseFilepath + "_" + time.toStringIO() + ".nc";
      io.read(coeffsSaver_.at(time), filepath);
    }
  }
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::read() done" << std::endl;
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::write() const {
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::write() starting" << std::endl;
  util::Timer timer(classname(), "write");
  eckit::LocalConfiguration outputConfig(coeffsConfig_, "output");
  const std::string baseFilepath = outputConfig.getString("base filepath");
  if (outputConfig.getBool("legacy", true)) {
    for (const auto & element : coeffsSaver_) {
      const std::string filepath = baseFilepath + "_" + element.first.toStringIO();
      outputConfig.set("filepath", filepath);
      util::writeFieldSet(updateGeometry_.getComm(), outputConfig, element.second);
    }
  } else {
    util::ParallelFieldSetIO io(updateGeometry_.functionSpace(),
                                outputConfig.getString("grid name"),
                                util::ParallelFieldSetIO::Mode::Write);
    for (const auto & element : coeffsSaver_) {
      const std::string filepath = baseFilepath + "_" + element.first.toStringIO() + ".nc";
      io.write(element.second, filepath);
    }
  }
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::write() done" << std::endl;
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::updateIncTL(Increment_ & dx) const {
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::updateIncTL() starting" << std::endl;
  const auto updateStencilArray = atlas::array::make_view<int, 2>(updateStencil_);
  atlas::FieldSet & dxFSet = dx.fieldSet().fieldSet();
  auto updateBuffer = util::perThreadStorage<double>(updateVars_.size()*nLevels_);
  auto bufferView = atlas::array::make_view<double, 2>(updateBuffer);
  util::IndexSpace1D owned_range = {0, static_cast<atlas::idx_t>(owned_.size())};
  util::for_each_index(owned_range,
  [&](atlas::idx_t idx) {
    const atlas::idx_t i = owned_[idx];
    auto updateValsView = bufferView.slice(atlas_omp_get_thread_num(),
      atlas::array::Range::all());
    updateValsView.assign(0.0);
    // Calculate update values
    for (size_t v = 0; v < updateVars_.size(); v++) {
      auto coeffsView  = atlas::array::make_view<double, 3>(
           coeffsSaver_.at(dx.validTime())[updateVars_[v].name()]);
      for (auto k = 0; k < nLevels_; k++) {
        for (size_t v2 = 0; v2 < updateVars_.size(); v2++) {
          auto dxArray = atlas::array::make_view<double, 2>(dxFSet[updateVars_[v2].name()]);
          for (auto s = 0; s < influenceSize_; s++) {
            updateValsView(k + v * nLevels_)
              += coeffsView(i, k, v2 * influenceSize_ + s) * dxArray(i, updateStencilArray(k, s));
          }
        }
      }
    }
    // Update column
    for (size_t v = 0; v < updateVars_.size(); v++) {
      auto dxArray = atlas::array::make_view<double, 2>(dxFSet[updateVars_[v].name()]);
      for (auto k = 0; k < nLevels_; k++) {
        dxArray(i, k) += updateValsView(k + v * nLevels_);
      }
    }
  });
  dx.synchronizeFields();
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::updateIncTL() done" << std::endl;
}

//------------------------------------------------------------------------------

template<typename MODEL>
void HybridLinearModelCoeffs<MODEL>::updateIncAD(Increment_ & dx) const {
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::updateIncAD() starting" << std::endl;
  const auto updateStencilArray = atlas::array::make_view<int, 2>(updateStencil_);
  atlas::FieldSet & dxFSet = dx.fieldSet().fieldSet();
  auto updateBuffer = util::perThreadStorage<double>(updateVars_.size()*nLevels_);
  auto bufferView = atlas::array::make_view<double, 2>(updateBuffer);
  util::IndexSpace1D owned_range = {0, static_cast<atlas::idx_t>(owned_.size())};
  util::for_each_index(owned_range,
  [&](atlas::idx_t idx) {
    const atlas::idx_t i = owned_[idx];
    auto updateValsView = bufferView.slice(atlas_omp_get_thread_num(),
      atlas::array::Range::all());
    updateValsView.assign(0.0);
    // Adjoint of "Update column"
    for (size_t v = 0; v < updateVars_.size(); v++) {
      auto dxArray = atlas::array::make_view<double, 2>(dxFSet[updateVars_[v].name()]);
      for (auto k = 0; k < nLevels_; k++) {
        updateValsView[k + v * nLevels_] += dxArray(i, k);
      }
    }
    // Adjoint of "Calculate update values"
    for (size_t v = 0; v < updateVars_.size(); v++) {
      auto coeffsView = atlas::array::make_view<double, 3>(
           coeffsSaver_.at(dx.validTime())[updateVars_[v].name()]);
      for (auto k = 0; k < nLevels_; k++) {
        for (size_t v2 = 0; v2 < updateVars_.size(); v2++) {
          auto dxArray = atlas::array::make_view<double, 2>(dxFSet[updateVars_[v2].name()]);
          for (auto s = 0; s < influenceSize_; s++) {
            dxArray(i, updateStencilArray(k, s))
              += coeffsView(i, k, v2 * influenceSize_ + s) * updateValsView(k + v * nLevels_);
          }
        }
      }
    }
  });
  dx.synchronizeFields();
  Log::trace() << "HybridLinearModelCoeffs<MODEL>::updateIncAD() done" << std::endl;
}

}  // namespace oops

#endif  // OOPS_GENERIC_HYBRIDLINEARMODELCOEFFS_H_
