/*
 * (C) Copyright 2023-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>

#include "oops/base/PostBase.h"
#include "oops/base/State.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

template <typename MODEL>
class JqTerm : public PostBase<State<MODEL>> {
  typedef State<MODEL>              State_;

 public:
  JqTerm();

  State_ & getMofX() {return *mx_;}

 private:
  void doProcessing(const State_ &) override {}
  void doFinalize(const State_ &) override;

// Data
  std::unique_ptr<State_> mx_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
JqTerm<MODEL>::JqTerm() : PostBase<State_>(), mx_() {}

// -----------------------------------------------------------------------------

template <typename MODEL>
void JqTerm<MODEL>::doFinalize(const State_ & mx) {
  Log::trace() << "JqTerm::doFinalize start" << std::endl;
  mx_ = std::make_unique<State_>(mx);
  Log::trace() << "JqTerm::doFinalize done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops
