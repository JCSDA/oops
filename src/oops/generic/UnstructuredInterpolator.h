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

#include "oops/generic/LocalInterpolatorBase.h"
#include "oops/util/ObjectCounter.h"

namespace oops {

class GeometryData;
class Variables;

template <typename MODEL>
class Geometry;
template <typename MODEL>
class Increment;
template <typename MODEL>
class State;

// -----------------------------------------------------------------------------

class UnstructuredInterpolator : public LocalInterpolatorBase,
                                 private util::ObjectCounter<UnstructuredInterpolator> {
 public:
  static const std::string classname() {return "oops::UnstructuredInterpolator";}

  UnstructuredInterpolator(const eckit::Configuration &, const GeometryData &,
                           const std::vector<double> &, const std::vector<double> &);

  // Constructor taking a MODEL-specific Geometry
  template <typename MODEL>
  UnstructuredInterpolator(const eckit::Configuration &, const Geometry<MODEL> &,
                           const std::vector<double> &, const std::vector<double> &);

  // Interpolator interface with no target-point mask, i.e., interpolates to every target point.
  void apply(const Variables &, const atlas::FieldSet &, std::vector<double> &) const override;
  void applyAD(const Variables &, atlas::FieldSet &, const std::vector<double> &) const override;

  // Interpolator interface with a target-point mask, i.e., interpolates to target points for which
  // the mask is true. At points for which mask is false, the return vector is unmodified from its
  // input state.
  void apply(const Variables &, const atlas::FieldSet &, const std::vector<bool> &,
             std::vector<double> &) const override;
  void applyAD(const Variables &, atlas::FieldSet &, const std::vector<bool> &,
               const std::vector<double> &) const override;

  // MODEL-specific interface to generic interpolator
  template <typename MODEL>
  void apply(const Variables &, const State<MODEL> &, std::vector<double> &) const;
  template <typename MODEL>
  void apply(const Variables &, const Increment<MODEL> &, std::vector<double> &) const;
  template <typename MODEL>
  void applyAD(const Variables &, Increment<MODEL> &, const std::vector<double> &) const;
  template <typename MODEL>
  void apply(const Variables &, const State<MODEL> &, const std::vector<bool> &,
             std::vector<double> &) const;
  template <typename MODEL>
  void apply(const Variables &, const Increment<MODEL> &, const std::vector<bool> &,
             std::vector<double> &) const;
  template <typename MODEL>
  void applyAD(const Variables &, Increment<MODEL> &, const std::vector<bool> &,
               const std::vector<double> &) const;

 private:
  // Small struct to help organize the interpolation matrices (= stencils and weights)
  struct InterpMatrix {
    std::vector<bool> targetHasValidStencil;
    std::vector<std::vector<size_t>> stencils;
    std::vector<std::vector<double>> weights;
  };

  void applyPerLevel(const InterpMatrix &,
                     const std::string &,
                     const std::vector<bool> &,
                     const atlas::array::ArrayView<double, 2> &,
                     std::vector<double>::iterator &, const size_t &) const;
  void applyPerLevelAD(const InterpMatrix &,
                       const std::string &,
                       const std::vector<bool> &,
                       atlas::array::ArrayView<double, 2> &,
                       std::vector<double>::const_iterator &, const size_t &) const;
  void print(std::ostream &) const override;

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

template <typename MODEL>
UnstructuredInterpolator::UnstructuredInterpolator(const eckit::Configuration & config,
                                                   const Geometry<MODEL> & geom,
                                                   const std::vector<double> & lats_out,
                                                   const std::vector<double> & lons_out)
  : UnstructuredInterpolator(config, geom.generic(), lats_out, lons_out) {}

// -----------------------------------------------------------------------------

template<typename MODEL>
void UnstructuredInterpolator::apply(const Variables & vars, const State<MODEL> & xx,
                                     std::vector<double> & locvals) const
{
  std::vector<bool> target_mask(nout_, true);
  this->apply(vars, xx.fieldSet().fieldSet(), target_mask, locvals);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void UnstructuredInterpolator::apply(const Variables & vars, const Increment<MODEL> & dx,
                                     std::vector<double> & locvals) const
{
  std::vector<bool> target_mask(nout_, true);
  this->apply(vars, dx.fieldSet().fieldSet(), target_mask, locvals);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void UnstructuredInterpolator::applyAD(const Variables & vars, Increment<MODEL> & dx,
                                       const std::vector<double> & vals) const {
  std::vector<bool> target_mask(nout_, true);
  this->applyAD(vars, dx.fieldSet().fieldSet(), target_mask, vals);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void UnstructuredInterpolator::apply(const Variables & vars, const State<MODEL> & xx,
                                     const std::vector<bool> & target_mask,
                                     std::vector<double> & locvals) const
{
  this->apply(vars, xx.fieldSet().fieldSet(), target_mask, locvals);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void UnstructuredInterpolator::apply(const Variables & vars, const Increment<MODEL> & dx,
                                     const std::vector<bool> & target_mask,
                                     std::vector<double> & locvals) const
{
  this->apply(vars, dx.fieldSet().fieldSet(), target_mask, locvals);
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void UnstructuredInterpolator::applyAD(const Variables & vars, Increment<MODEL> & dx,
                                       const std::vector<bool> & target_mask,
                                       const std::vector<double> & vals) const {
  this->applyAD(vars, dx.fieldSet().fieldSet(), target_mask, vals);
}

// -----------------------------------------------------------------------------

}  // namespace oops
