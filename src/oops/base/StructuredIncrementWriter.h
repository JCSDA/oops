/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/PostBase.h"
#include "oops/base/State.h"
#include "oops/base/StructuredGridWriter.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class StructuredIncrementWriter : public PostBase<State<MODEL> > {
  typedef Geometry<MODEL>    Geometry_;
  typedef Increment<MODEL>   Increment_;
  typedef State<MODEL>       State_;
  typedef StructuredGridWriter<MODEL> StructuredGridWriter_;

 public:
  explicit StructuredIncrementWriter(const eckit::Configuration & conf);
  ~StructuredIncrementWriter() = default;

 private:
  std::unique_ptr<Geometry_> resol_;
  const eckit::LocalConfiguration conf_;
  int pass_;
  std::vector<std::unique_ptr<State_>> bgs_;
  typename std::vector<std::unique_ptr<State_>>::iterator current_;
  std::unique_ptr<StructuredGridWriter_> writer_;

  void doInitialize(const State_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessing(const State_ &) override;
  void doFinalize(const State_ &) override;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
StructuredIncrementWriter<MODEL>::StructuredIncrementWriter(const eckit::Configuration & conf)
  : PostBase<State_>(conf), resol_(), conf_(conf), pass_(0) {}

// -----------------------------------------------------------------------------

template <typename MODEL>
void StructuredIncrementWriter<MODEL>::doInitialize(const State_ & xx,
                                          const util::DateTime &, const util::Duration &) {
  pass_ += 1;
  ASSERT(pass_ == 1 || pass_ == 2);
  if (pass_ == 1) {
    eckit::LocalConfiguration resConf(conf_, "geometry");
    resol_ = std::make_unique<Geometry_>(resConf, xx.geometry().getComm(),
                                                  xx.geometry().timeComm());
  }
  if (pass_ == 2) {
    writer_ = std::make_unique<StructuredGridWriter_>(conf_, *resol_);
    current_ = bgs_.begin();
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void StructuredIncrementWriter<MODEL>::doProcessing(const State_ & xx) {
  if (pass_ == 1) {
    bgs_.emplace_back(std::make_unique<State_>(*resol_, xx));
  } else {
    ASSERT(xx.validTime() == (*current_)->validTime());
    ASSERT(*current_);
    State_ xlr(*resol_, xx);
    Increment_ dx(*resol_, xx.variables(), xx.validTime());
    dx.diff(xlr, **current_);
    Log::info() << "Writing Increment to Structured grid:" << dx << std::endl;
    Log::test() << "Writing Increment to Structured grid:" << dx << std::endl;
    writer_->interpolateAndWrite(dx, xlr);
    current_->reset();
    ++current_;
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void StructuredIncrementWriter<MODEL>::doFinalize(const State_ &) {
  if (pass_ == 2) {
    bgs_.clear();
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops
