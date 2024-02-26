/*
 * (C) Copyright 2020-2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_INTERFACE_ANALYTICINITBASE_H_
#define OOPS_INTERFACE_ANALYTICINITBASE_H_

#include "oops/generic/AnalyticInitBase.h"

namespace oops {

namespace interface {

/// \brief Initializes GeoVaLs with analytic formula
template <typename OBS>
class AnalyticInitBase : public oops::AnalyticInitBase<OBS> {
  typedef typename OBS::GeoVaLs          GeoVaLs_;
  typedef typename OBS::SampledLocations SampledLocations_;

 public:
  // Overrides of oops::AnalyticInit methods, converting between oops interface classes and
  // OBS-specific classes operated upon by OBS-specific implementations of the AnalyticInit
  // interface.
  void fillGeoVaLs(const SampledLocations<OBS> & locs, GeoVaLs<OBS> & geovals) const final {
    fillGeoVaLs(locs.sampledLocations(), geovals.geovals());
  }

/*! \brief GeoVaLs Analytic Initialization
 *
 * \details **AnalyticInit::fillGeoVaLs** was introduced in May, 2018 (initially
 * as a GeoVaLs constructor) for use with the interpolation test. The interpolation test
 * requires an initialization of a GeoVaLs object based on the same analytic
 * formulae used for the State initialization (see test::TestStateInterpolation()
 * for further information).  This in turn requires information about the
 * vertical profile in addition to the latitude and longitude positional
 * information in the SampledLocations object.  Currently, this information
 * about the vertical profile is obtained from an existing GeoVaLs object
 * (passed as *gvals*) that represents the output of the State::interpolate()
 * method.
 *
 * \date May, 2018: created as a constructor (M. Miesch, JCSDA)
 * \date June, 2018: moved to a method (M. Miesch, JCSDA)
 *
 * \sa test::TestStateInterpolation()
 */
  virtual void fillGeoVaLs(const SampledLocations_ &, GeoVaLs_ &) const = 0;
};

}  // namespace interface

}  // namespace oops

#endif  // OOPS_INTERFACE_ANALYTICINITBASE_H_
