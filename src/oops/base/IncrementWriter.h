/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <vector>

#include "eckit/config/Configuration.h"

#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/PostBase.h"
#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class IncrementWriter : public PostBase<State<MODEL> > {
  typedef Geometry<MODEL>    Geometry_;
  typedef Increment<MODEL>   Increment_;
  typedef State<MODEL>       State_;

 public:
  explicit IncrementWriter(const eckit::Configuration & conf);
  ~IncrementWriter() = default;

 private:
  std::unique_ptr<Geometry_> resol_;
  const Variables vars_;
  const eckit::LocalConfiguration conf_;
  int pass_;
  std::vector<std::unique_ptr<State_>> bgs_;
  typename std::vector<std::unique_ptr<State_>>::iterator current_;

  void doInitialize(const State_ &, const util::DateTime &, const util::Duration &) override;
  void doProcessing(const State_ &) override;
  void doFinalize(const State_ &) override;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
IncrementWriter<MODEL>::IncrementWriter(const eckit::Configuration & conf)
  : PostBase<State_>(conf), resol_(), vars_(conf, "variables"), conf_(conf), pass_(0) {}

// -----------------------------------------------------------------------------

template <typename MODEL>
void IncrementWriter<MODEL>::doInitialize(const State_ & xx,
                                          const util::DateTime &, const util::Duration &) {
  pass_ += 1;
  ASSERT(pass_ == 1 || pass_ == 2);
  if (pass_ == 1 && conf_.has("geometry")) {
    eckit::LocalConfiguration resConf(conf_, "geometry");
    resol_ = std::make_unique<Geometry_>(resConf, xx.geometry().getComm(),
                                                  xx.geometry().timeComm());
  }
  if (pass_ == 2) current_ = bgs_.begin();
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void IncrementWriter<MODEL>::doProcessing(const State_ & xx) {
  if (pass_ == 1) {
    if (resol_) {
      bgs_.emplace_back(std::make_unique<State_>(*resol_, xx));
    } else {
      bgs_.emplace_back(std::make_unique<State_>(xx));
    }
  } else {
    ASSERT(xx.validTime() == (*current_)->validTime());
    ASSERT(*current_);
    if (resol_) {
      State_ xlr(*resol_, xx);
      Increment_ dx(*resol_, vars_, xx.validTime());
      dx.diff(xlr, **current_);
      Log::info() << "Writing Increment:" << dx << std::endl;
      Log::test() << "Writing Increment:" << dx << std::endl;
      dx.write(conf_);
      current_->reset();
    } else {
      Increment_ dx(xx.geometry(), vars_, xx.validTime());
      dx.diff(xx, **current_);
      Log::info() << "Writing Increment:" << dx << std::endl;
      Log::test() << "Writing Increment:" << dx << std::endl;
      dx.write(conf_);
      current_->reset();
    }
    ++current_;
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void IncrementWriter<MODEL>::doFinalize(const State_ &) {
  if (pass_ == 2) {
    bgs_.clear();
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops
