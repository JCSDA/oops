/*
 * (C) Copyright 2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "oops/base/PostBase.h"

#include "eckit/config/LocalConfiguration.h"
#include "oops/util/DateTime.h"
#include "oops/util/workflow.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename FLDS> class WorkflowUpdater : public PostBase<FLDS> {
 public:
  explicit WorkflowUpdater(const eckit::Configuration & conf): PostBase<FLDS>(conf) {}
  ~WorkflowUpdater() {}

 private:
  void doProcessing(const FLDS & xx) override;
  void doInitialize(const FLDS &, const util::DateTime &, const util::Duration &) override;

  util::DateTime start_;
};

// -----------------------------------------------------------------------------

template<typename FLDS>
void WorkflowUpdater<FLDS>::doInitialize(const FLDS & xx,
                                         const util::DateTime &, const util::Duration &) {
  start_ = xx.validTime();
}

// -----------------------------------------------------------------------------

template<typename FLDS>
void WorkflowUpdater<FLDS>::doProcessing(const FLDS & xx) {
  const util::Duration current = xx.validTime() - start_;
  const int hh = round(current.toSeconds() / 3600);
  Log::info() << "WorkflowUpdater: step = " << current << " " << hh  << std::endl;
  util::update_workflow_meter("step", hh);
}

// -----------------------------------------------------------------------------

}  // namespace oops

