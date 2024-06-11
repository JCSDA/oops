/*
 * (C) Crown copyright 2023, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_LOCATIONS_H_
#define OOPS_BASE_LOCATIONS_H_

#include <memory>
#include <utility>
#include <vector>

#include "oops/base/SamplingMethodSelector.h"
#include "oops/base/TrivialSamplingMethodSelector.h"
#include "oops/base/Variable.h"
#include "oops/interface/SampledLocations.h"

namespace oops {

/// \brief Observation locations.
///
/// OOPS uses observation locations when computing GeoVaLs (model variables interpolated at
/// observation locations).
///
/// In general, observation locations, i.e. the regions of space probed by individual observations,
/// may have complex geometry. To compute GeoVaLs or ObsDiagnostics, OOPS does not need to know
/// that geometry exactly. Instead, the observation locations need to be discretized by sampling
/// them with a set of interpolation paths (1D rays along which the model fields will be
/// interpolated). Currently these paths must be vertical, but there are plans to add support for
/// slanted paths in the future. Different sampling methods may be used for different variables:
/// for instance, for rapidly varying variables each location may need to be sampled with multiple
/// interpolation paths, whereas for variables varying more slowly a single path per location may
/// be enough. The set of paths produced by each sampling method is represented by a
/// SampledLocations object. The Locations class stores one or more such objects and maps each
/// variable to one of them.
template <typename OBS>
class Locations {
 public:
  typedef SampledLocations<OBS> SampledLocations_;

  /// \brief Creates an object holding a single set of interpolation paths sampling the observation
  /// locations, to be used for all variables.
  ///
  /// \note Keeping this constructor implicit simplifies code in the common case where a single
  /// discretization of the observation locations is enough.
  // NOLINTNEXTLINE(runtime/explicit)
  Locations(SampledLocations_ sampledLocations) {
    sampledLocations_.push_back(std::move(sampledLocations));
    selector_ = std::make_unique<TrivialSamplingMethodSelector>();
  }

  /// \brief Creates an object holding one or more sets of interpolation paths sampling the
  /// observation locations. `selector` is used to map each variable to one of these sets.
  Locations(std::vector<SampledLocations_> sampledLocations,
            std::unique_ptr<SamplingMethodSelector> selector)
    : sampledLocations_(std::move(sampledLocations)), selector_(std::move(selector))
  {}

  /// \brief Returns the number of methods by which the observation locations have been sampled.
  size_t numSamplingMethods() const { return sampledLocations_.size(); }

  /// \brief Returns the set of interpolation paths produced by the `i`th location sampling method.
  const SampledLocations_& samplingMethod(size_t i) const { return sampledLocations_.at(i); }

  /// \brief Returns the set of paths along which the variable `varName` should be interpolated.
  const SampledLocations_& samplingMethod(const Variable &varName) const {
    return sampledLocations_[samplingMethodIndex(varName)];
  }

  /// \brief Returns the index of the location sampling method defining the set of paths along
  /// which the variable `varName` should be interpolated.
  size_t samplingMethodIndex(const Variable &varName) const {
    const size_t result = selector_->methodIndex(varName);
    ASSERT(result < sampledLocations_.size());
    return result;
  }

 private:
  std::vector<SampledLocations_> sampledLocations_;
  std::unique_ptr<SamplingMethodSelector> selector_;
};

}  // namespace oops

#endif  // OOPS_BASE_LOCATIONS_H_
