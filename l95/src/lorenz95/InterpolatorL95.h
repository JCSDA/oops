/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <vector>

#include "oops/util/Printable.h"

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
}

namespace lorenz95 {
  class Resolution;
  class IncrementL95;
  class StateL95;

// -----------------------------------------------------------------------------

class InterpolatorL95 : public util::Printable {
 public:
  InterpolatorL95(const eckit::Configuration &, const Resolution &,
                  const std::vector<double> &, const std::vector<double> &);
  ~InterpolatorL95();

  void apply(const oops::Variables &, const StateL95 &, const std::vector<bool> &,
             std::vector<double> &) const;
  void apply(const oops::Variables &, const IncrementL95 &, const std::vector<bool> &,
             std::vector<double> &) const;
  void applyAD(const oops::Variables &, IncrementL95 &, const std::vector<bool> &,
               const std::vector<double> &) const;

 private:
  void print(std::ostream &) const;

  size_t nout_;
  std::vector<size_t> ilocs_;
};
// -----------------------------------------------------------------------------

}  // namespace lorenz95

