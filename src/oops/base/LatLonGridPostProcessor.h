/*
 * (C) Copyright 2022- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "eckit/config/Configuration.h"

#include "oops/base/Geometry.h"
#include "oops/base/LatLonGridWriter.h"
#include "oops/base/PostBase.h"

namespace oops {

// -----------------------------------------------------------------------------

/// PostProcessor that wraps the LatLonGridWriter
template <typename MODEL, typename FLDS>
class LatLonGridPostProcessor : public PostBase<FLDS>, public LatLonGridWriter<MODEL> {
 public:
  explicit LatLonGridPostProcessor(
      const eckit::Configuration & conf,
      const Geometry<MODEL> & sourceGeometry);
  ~LatLonGridPostProcessor() = default;

 private:
  void doProcessing(const FLDS & xx) override;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename FLDS>
LatLonGridPostProcessor<MODEL, FLDS>::LatLonGridPostProcessor(
    const eckit::Configuration & conf,
    const Geometry<MODEL> & sourceGeometry)
: PostBase<FLDS>(conf),
  LatLonGridWriter<MODEL>(conf, sourceGeometry) {}

// -----------------------------------------------------------------------------

template <typename MODEL, typename FLDS>
void LatLonGridPostProcessor<MODEL, FLDS>::doProcessing(const FLDS & xx) {
  // supports output on pressure levels for FLDS of State type only
  this->interpolateAndWrite(xx);
}

// -----------------------------------------------------------------------------

}  // namespace oops
