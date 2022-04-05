/*
 * (C) Copyright 2020-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "eckit/config/Configuration.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace oops {
  class Variables;

// -----------------------------------------------------------------------------

/// \brief Interface for Unstructured interpolation

class InterpolatorUnstructured : public util::Printable,
                                 private util::ObjectCounter<InterpolatorUnstructured> {
 public:
  static const std::string classname() {return "oops::InterpolatorUnstructured";}

  InterpolatorUnstructured(const eckit::Configuration &, const std::vector<double> &,
                           const std::vector<double> &, const std::vector<double> &);
  ~InterpolatorUnstructured() {}

  void apply(const Variables &, const atlas::FieldSet &, std::vector<double> &) const;
  void applyAD(const Variables &, atlas::FieldSet &, const std::vector<double> &) const;

 private:
  std::string interp_method_;
  int nninterp_;
  size_t nout_;
  std::vector<std::vector<size_t>> interp_i_;
  std::vector<std::vector<double>> interp_w_;

  void apply1lev(const std::string &, const atlas::array::ArrayView<double, 1> &,
                 std::vector<double>::iterator &) const;
  void applyLevs(const std::string &, const atlas::array::ArrayView<double, 2> &,
                 std::vector<double>::iterator &, const size_t &) const;
  void apply1levAD(const std::string &, atlas::array::ArrayView<double, 1> &,
                   std::vector<double>::const_iterator &) const;
  void applyLevsAD(const std::string &, atlas::array::ArrayView<double, 2> &,
                   std::vector<double>::const_iterator &, const size_t &) const;
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace oops
