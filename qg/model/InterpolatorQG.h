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

namespace qg {
  class GeometryQG;
  class IncrementQG;
  class StateQG;

// -----------------------------------------------------------------------------

class InterpolatorQG : public util::Printable {
 public:
  InterpolatorQG(const eckit::Configuration &, const GeometryQG &, const std::vector<double> &);
  ~InterpolatorQG();

  void apply(const oops::Variables &, const StateQG &, std::vector<double> &) const;
  void apply(const oops::Variables &, const IncrementQG &, std::vector<double> &) const;
  void applyAD(const oops::Variables &, IncrementQG &, const std::vector<double> &) const;

 private:
  void print(std::ostream &) const;

  const size_t nlevs_;
  const size_t nlocs_;
  std::vector<double> locs_;
};
// -----------------------------------------------------------------------------

}  // namespace qg

