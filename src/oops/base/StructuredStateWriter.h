/*
 * (C) Copyright 2022-2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <string>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/PostBase.h"
#include "oops/base/State.h"
#include "oops/base/StructuredGridWriter.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class StructuredStateWriter : public PostBase<State<MODEL>> {
  typedef State<MODEL>     State_;
  typedef StructuredGridWriter<MODEL> StructuredGridWriter_;
 public:
  explicit StructuredStateWriter(const eckit::Configuration & conf);
  ~StructuredStateWriter() = default;

 private:
  eckit::LocalConfiguration conf_;
  std::unique_ptr<StructuredGridWriter_> writer_;

  void doInitialize(const State_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessing(const State_ &) override;
  void doFinalize(const State_ &) override;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
StructuredStateWriter<MODEL>::StructuredStateWriter(const eckit::Configuration & conf)
: PostBase<State_>(conf), conf_(conf) {}

// -----------------------------------------------------------------------------

template <typename MODEL>
void StructuredStateWriter<MODEL>::doInitialize(const State_ & xx, const util::DateTime &,
                                                const util::Duration &) {
  writer_ = std::make_unique<StructuredGridWriter_>(conf_, xx.geometry());
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void StructuredStateWriter<MODEL>::doProcessing(const State_ & xx) {
  writer_->interpolateAndWrite(xx);
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void StructuredStateWriter<MODEL>::doFinalize(const State_ &) {
  writer_.reset();
}

// -----------------------------------------------------------------------------

}  // namespace oops
