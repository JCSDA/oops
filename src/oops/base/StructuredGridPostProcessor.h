/*
 * (C) Copyright 2022-2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "eckit/config/Configuration.h"

#include "oops/base/Geometry.h"
#include "oops/base/PostBase.h"
#include "oops/base/StructuredGridWriter.h"

namespace oops {

// -----------------------------------------------------------------------------

/// PostProcessor that wraps the StructuredGridWriter
template <typename MODEL, typename FLDS>
class StructuredGridPostProcessor : public PostBase<FLDS>, public StructuredGridWriter<MODEL> {
 public:
  explicit StructuredGridPostProcessor(
      const eckit::Configuration & conf,
      const Geometry<MODEL> & sourceGeometry);
  ~StructuredGridPostProcessor() = default;

 private:
  void doProcessing(const FLDS & xx) override;
};

// -----------------------------------------------------------------------------

template <typename MODEL, typename FLDS>
StructuredGridPostProcessor<MODEL, FLDS>::StructuredGridPostProcessor(
    const eckit::Configuration & conf,
    const Geometry<MODEL> & sourceGeometry)
: PostBase<FLDS>(conf),
  StructuredGridWriter<MODEL>(conf, sourceGeometry) {}

// -----------------------------------------------------------------------------

template <typename MODEL, typename FLDS>
void StructuredGridPostProcessor<MODEL, FLDS>::doProcessing(const FLDS & xx) {
  // supports output on pressure levels for FLDS of State type only
  this->interpolateAndWrite(xx);
}

// -----------------------------------------------------------------------------

}  // namespace oops
