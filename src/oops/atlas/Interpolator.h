/*
 * (C) Copyright 2023-2025 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <vector>

#include "atlas/field.h"
#include "eckit/config/Configuration.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/generic/SourceProximityPartitioner.h"


namespace oops {
namespace atlasbase {

// Base class describing an interpolator that satisfies the LocalInterpolator
// interface but is implemented to act on Atlas data structures.
class Interpolator {
 public:
  Interpolator(const eckit::Configuration & conf,
               const GeometryData & source_geometry,
               const std::vector<double> & target_lats,
               const std::vector<double> & target_lons);
  virtual ~Interpolator() = default;

  // Partition target points based on proximity to a source grid node.
  static SourceProximityPartitioner makeTargetPartitioner(const GeometryData & geom) {
    return makeSourceProximityPartitioner(geom.functionSpace(), geom.comm());
  }

  // Apply a halo exchange to the atlas::FieldSet before interpolating.
  // This preprocess action is associated with the use of atlas::FieldSets to
  // represent the model fields, so shouldn't depend on the interpolator
  // specialization or internal state.
  static void preprocess(atlas::FieldSet &);
  static void preprocessAD(atlas::FieldSet &);

  // Implement an Atlas-based interpolation algorithm; must be provided by each
  // derived class.
  //
  // Note,
  // - `mask` is a vector over the target points; the interpolation should only
  //   be computed for target points whose corresponding mask element is true.
  // - `buffer` is a contiguous vector into which the interpolation results are
  //   written to, typically for redistribution by the client code.
  //   - The buffer should be resized to (number of target points * number of
  //     field levels), where `field levels` is the total number of levels
  //     across all variables, keeping in mind that different variables in the
  //     FieldSet may have different numbers of levels.
  //   - The buffer ordering should be:
  //     - fastest-varying (contiguous) index is level
  //     - middle index is target point
  //     - slowest-varying index is variable
  //     The buffer can be thought of as columns of data appended together,
  //     matching the data ordering of the atlas::FieldSet and the UFO GeoVaLs.
  //   - When a target point is excluded by the mask (mask element == false),
  //     the corresponding buffer entries should not be written to. This permits
  //     successive calls to the interpolation, using different masks, to fill
  //     in different portions of the buffer.
  virtual void apply(const Variables& vars,
                     const atlas::FieldSet& fields,
                     const std::vector<bool>& mask,
                     std::vector<double>& buffer) const = 0;
  virtual void applyAD(const Variables& vars,
                       atlas::FieldSet& fields,
                       const std::vector<bool>& mask,
                       const std::vector<double>& buffer) const = 0;

  // Convenience wrappers that provide an unmasked interpolation interface.
  void apply(const Variables& vars,
             const atlas::FieldSet& fields,
             std::vector<double>& buffer) const;
  void applyAD(const Variables& vars,
               atlas::FieldSet& fields,
               const std::vector<double>& buffer) const;

  // Transfer the interpolation's output buffer back into a FieldSet.
  // These methods do NOT rely on any internal state of the interpolator, they
  // only encode the inverse of the transformation done in apply() to get an
  // output buffer from the FieldSet
  static void bufferToFieldSet(const Variables &,
                               const std::vector<size_t> &,
                               const std::vector<double> &,
                               atlas::FieldSet &);
  static void bufferToFieldSetAD(const Variables &,
                                 const std::vector<size_t> &,
                                 std::vector<double> &,
                                 const atlas::FieldSet &);

 private:
  int nb_targets_ = 0;
};

}  // namespace atlasbase
}  // namespace oops
