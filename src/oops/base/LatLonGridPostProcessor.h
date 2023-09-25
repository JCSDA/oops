/*
 * (C) Copyright 2022- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "oops/base/Geometry.h"
#include "oops/base/LatLonGridWriter.h"
#include "oops/base/PostBase.h"
#include "oops/util/parameters/Parameters.h"

namespace oops {

// -----------------------------------------------------------------------------

class LatLonGridPostProcessorParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(LatLonGridPostProcessorParameters, Parameters)

 public:
  // time steps at which the lat-lon grid data is written out
  PostTimerParameters postTimer{this};
  // output grid, fields, and file options
  LatLonGridWriterParameters latlonWriter{this};
};

// -----------------------------------------------------------------------------

/// PostProcessor that wraps the LatLonGridWriter
template <typename MODEL, typename FLDS>
class LatLonGridPostProcessor : public PostBase<FLDS>, public LatLonGridWriter<MODEL> {
 public:
  explicit LatLonGridPostProcessor(
      const LatLonGridPostProcessorParameters & parameters,
      const Geometry<MODEL> & sourceGeometry);
  ~LatLonGridPostProcessor() = default;

 private:
  void doProcessing(const FLDS & xx) override;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename FLDS>
LatLonGridPostProcessor<MODEL, FLDS>::LatLonGridPostProcessor(
    const LatLonGridPostProcessorParameters & parameters,
    const Geometry<MODEL> & sourceGeometry)
: PostBase<FLDS>(parameters.postTimer),
  LatLonGridWriter<MODEL>(parameters.latlonWriter, sourceGeometry) {}

// -----------------------------------------------------------------------------

template <typename MODEL, typename FLDS>
void LatLonGridPostProcessor<MODEL, FLDS>::doProcessing(const FLDS & xx) {
  this->interpolateAndWrite(xx);  // supports output on pressure levels for FLDS of State type only
}

// -----------------------------------------------------------------------------

}  // namespace oops
