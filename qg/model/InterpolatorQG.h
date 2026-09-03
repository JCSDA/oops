/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <vector>

#include "oops/generic/UnstructuredInterpolator.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Printable.h"

namespace atlas {
  class FieldSet;
}

namespace eckit {
  class Configuration;
}

namespace oops {
  class Variables;
  class SourceProximityPartitioner;
}

namespace qg {
  class FieldsQG;
  class GeometryQG;
  class IncrementQG;
  class StateQG;

// -----------------------------------------------------------------------------

class InterpolatorQG : public util::Printable {
 public:
  InterpolatorQG(const eckit::Configuration &, const GeometryQG &,
                 const std::vector<double> &, const std::vector<double> &);
  ~InterpolatorQG();

  // Partition target points based on proximity to a source grid point.
  static oops::SourceProximityPartitioner makeTargetPartitioner(const GeometryQG &);

  static void preprocess(StateQG &);
  static void preprocess(IncrementQG &);
  static void preprocessAD(IncrementQG &);

  void apply(const oops::Variables &, const StateQG &, const std::vector<bool> &,
             std::vector<double> &) const;
  void apply(const oops::Variables &, const IncrementQG &, const std::vector<bool> &,
             std::vector<double> &) const;
  void applyAD(const oops::Variables &, IncrementQG &, const std::vector<bool> &,
               const std::vector<double> &) const;

 private:
  void apply(const oops::Variables &, const FieldsQG &, const std::vector<bool> &,
             std::vector<double> &) const;
  void print(std::ostream &) const override;

  const GeometryQG & geom_;
  const size_t nlocs_;
  std::vector<double> locs_;
};
// -----------------------------------------------------------------------------

}  // namespace qg

