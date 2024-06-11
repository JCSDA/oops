/*
 * (C) Crown copyright 2023, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_TRIVIALSAMPLINGMETHODSELECTOR_H_
#define OOPS_BASE_TRIVIALSAMPLINGMETHODSELECTOR_H_

#include "oops/base/SamplingMethodSelector.h"

namespace oops {

class Variable;

/// \brief Selects sampling method #0 for all variables.
class TrivialSamplingMethodSelector : public SamplingMethodSelector {
 public:
  size_t methodIndex(const Variable &/*varName*/) const override { return 0; }
};

}  // namespace oops

#endif  // OOPS_BASE_TRIVIALSAMPLINGMETHODSELECTOR_H_
