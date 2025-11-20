/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>
#include <unordered_map>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "eckit/config/Configuration.h"

#include "oops/atlas/Interpolator.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace oops {

class GeometryData;
class Variables;

// -----------------------------------------------------------------------------

class UnstructuredInterpolator : public atlasbase::Interpolator,
                                 public util::Printable,
                                 private util::ObjectCounter<UnstructuredInterpolator> {
 public:
  static const std::string classname() {return "oops::UnstructuredInterpolator";}

  UnstructuredInterpolator(const eckit::Configuration &, const GeometryData &,
                           const std::vector<double> &, const std::vector<double> &);
  UnstructuredInterpolator() = delete;

  using atlasbase::Interpolator::apply;
  using atlasbase::Interpolator::applyAD;

  // Interpolator interface with a target-point mask, i.e., interpolates to target points for which
  // the mask is true. At points for which mask is false, the return vector is unmodified from its
  // input state.
  void apply(const Variables &, const atlas::FieldSet &, const std::vector<bool> &,
             std::vector<double> &) const override;
  void applyAD(const Variables &, atlas::FieldSet &, const std::vector<bool> &,
               const std::vector<double> &) const override;

 private:
  // Small struct to help organize the interpolation matrices (= stencils and weights)
  struct InterpMatrix {
    std::vector<bool> targetHasValidStencil;
    std::vector<std::vector<size_t>> stencils;
    std::vector<std::vector<double>> weights;
  };

  void print(std::ostream &) const override;

  void doApply(const InterpMatrix &,
               const std::string &,
               const std::vector<bool> &,
               const atlas::array::ArrayView<double, 2> &,
               atlas::array::ArrayView<double, 2> &) const;
  void doApplyAD(const InterpMatrix &,
                 const std::string &,
                 const std::vector<bool> &,
                 atlas::array::ArrayView<double, 2> &,
                 const atlas::array::ArrayView<double, 2> &) const;

  void computeUnmaskedInterpMatrix(std::vector<double>, std::vector<double>) const;
  void computeMaskedInterpMatrix(const std::string &,
                                 const atlas::array::ArrayView<double, 2> &) const;

  const GeometryData & geom_;
  size_t nout_;

  // Current triangulation-based algorithm requires the interpolation stencil to contain 3 points,
  // but this could in principle depend on the config if we generalize the algorithm to perform
  // higher-order (than linear) interpolations using information from neighboring triangles.
  constexpr static size_t nstencil_ = 3;

  // The interpolation matrices depend on the mask used at runtime. We cache the matrices as they
  // are computed, to save computations across multiple interpolations using the same mask.
  // The caching is an implementation detail, so is done using a mutable member to preserve a
  // const interpolation interface. This may break threadsafety!
  mutable std::unordered_map<std::string, InterpMatrix> interp_matrices_;
  const std::string unmaskedName_{"unmasked"};
};

// -----------------------------------------------------------------------------

}  // namespace oops
