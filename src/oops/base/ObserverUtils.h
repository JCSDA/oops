/*
 * (C) Crown Copyright 2022-2023 Met Office UK.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */
#ifndef OOPS_BASE_OBSERVERUTILS_H_
#define OOPS_BASE_OBSERVERUTILS_H_

#include <algorithm>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "oops/base/Geometry.h"
#include "oops/base/Locations.h"
#include "oops/base/Variables.h"
#include "oops/interface/GeoVaLs.h"

namespace oops {

template <typename MODEL, typename OBS> class GetValues;

/// \brief Identify groups of variables interpolated along paths produced by the same observation
/// location sampling method.
template <typename OBS>
std::vector<Variables> groupVariablesByLocationSamplingMethod(
    const Variables &vars,
    const Locations<OBS> &locations) {
  std::vector<Variables> groupedVars(locations.numSamplingMethods());
  for (size_t v = 0; v < vars.size(); ++v)
    groupedVars[locations.samplingMethodIndex(vars[v])].push_back(vars[v]);
  return groupedVars;
}

/// \brief Return the number of levels in each GeoVaL
template<typename MODEL>
std::vector<std::vector<size_t>> variableSizes(const std::vector<Variables> &vars,
                                               const Geometry<MODEL> &geom) {
  std::vector<std::vector<size_t>> sizes;
  sizes.reserve(vars.size());
  std::transform(vars.begin(), vars.end(), std::back_inserter(sizes),
                 [&geom](const Variables &vs) { return geom.variableSizes(vs); });
  return sizes;
}

inline std::pair<Variables, std::vector<size_t>> mergeVariablesAndSizes(
    const std::vector<Variables> &groupedVars,
    const std::vector<std::vector<size_t>> &groupedSizes) {
  Variables mergedVars;
  // We use push_back rather than += to preserve order
  for (const Variables &vars : groupedVars)
    for (size_t i = 0; i != vars.size(); ++i)
      mergedVars.push_back(vars[i]);

  std::vector<size_t> mergedSizes;
  for (const std::vector<size_t> &sizes : groupedSizes)
    for (size_t i = 0; i != sizes.size(); ++i)
      mergedSizes.push_back(sizes[i]);

  return std::make_pair(std::move(mergedVars), std::move(mergedSizes));
}

template <typename MODEL, typename OBS>
std::vector<std::shared_ptr<GetValues<MODEL, OBS>>> makeGetValuesVector(
    const eckit::Configuration & conf, const Geometry<MODEL> & geom,
    const util::TimeWindow & timeWindow,
    const Locations<OBS> & locations,
    const std::vector<Variables> & varsGroupedBySamplingMethod) {
  typedef GetValues<MODEL, OBS> GetValues_;
  ASSERT(locations.numSamplingMethods() == varsGroupedBySamplingMethod.size());

  std::vector<std::shared_ptr<GetValues_>> getvalues;
  getvalues.reserve(locations.numSamplingMethods());
  for (size_t m = 0; m < locations.numSamplingMethods(); ++m) {
    getvalues.push_back(std::make_shared<GetValues_>(conf, geom, timeWindow,
                                                     locations.samplingMethod(m),
                                                     varsGroupedBySamplingMethod[m]));
  }
  return getvalues;
}

template <typename MODEL, typename OBS>
std::vector<std::shared_ptr<GetValues<MODEL, OBS>>> makeGetValuesVector(
    const eckit::Configuration & conf, const Geometry<MODEL> & geom,
    const util::TimeWindow & timeWindow,
    const Locations<OBS> & locations,
    const std::vector<Variables> & hopVarsGroupedBySamplingMethod,
    const std::vector<Variables> & hoptladVarsGroupedBySamplingMethod) {

  typedef GetValues<MODEL, OBS> GetValues_;
  ASSERT(locations.numSamplingMethods() == hopVarsGroupedBySamplingMethod.size());
  ASSERT(locations.numSamplingMethods() == hoptladVarsGroupedBySamplingMethod.size());

  std::vector<std::shared_ptr<GetValues_>> getvalues;
  getvalues.reserve(locations.numSamplingMethods());
  for (size_t m = 0; m < locations.numSamplingMethods(); ++m) {
    getvalues.push_back(std::make_shared<GetValues_>(
                          conf, geom, timeWindow, locations.samplingMethod(m),
                          hopVarsGroupedBySamplingMethod[m],
                          hoptladVarsGroupedBySamplingMethod[m]));
  }
  return getvalues;
}

template <typename MODEL, typename OBS>
GeoVaLs<OBS> makeAndFillGeoVaLs(
    const Locations<OBS> & locations,
    const Variables & vars, const std::vector<size_t> & varsizes,
    std::vector<std::shared_ptr<GetValues<MODEL, OBS>>> getvals) {

  GeoVaLs<OBS> geovals(locations, vars, varsizes);
  for (size_t m = 0; m < getvals.size(); ++m) {
    if (getvals[m]->useMethodsTL()) {
      getvals[m]->fillGeoVaLsTL(geovals);
    } else {
      getvals[m]->fillGeoVaLs(geovals);
    }
  }

  return geovals;
}

}  // namespace oops

#endif  // OOPS_BASE_OBSERVERUTILS_H_
