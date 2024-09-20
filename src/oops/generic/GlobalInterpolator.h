/*
 * (C) Copyright 2022-2023 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/Configuration.h"

#include "oops/generic/LocalInterpolatorBase.h"
#include "oops/mpi/mpi.h"
#include "oops/util/Printable.h"


namespace oops {

class GeometryData;

// -----------------------------------------------------------------------------

class GlobalInterpolator : public util::Printable {
 public:
  GlobalInterpolator(const eckit::Configuration &,
                     const GeometryData &,
                     const atlas::FunctionSpace &,
                     const eckit::mpi::Comm &);
  ~GlobalInterpolator() = default;

  /// Interpolate from source to target.
  ///
  /// The target FieldSet argument must either,
  /// - be empty, in which case it will be allocated within apply, OR
  /// - be correctly allocated, i.e., should contain one Field for each Field in source, with each
  ///   target Field allocated to the size of the target FunctionSpace and the levels of the
  ///   corresponding source Field.
  void apply(const atlas::FieldSet & source, atlas::FieldSet & target) const;

  /// Perform the adjoint of interpolation from source to target, moving data from target to source.
  ///
  /// The source FieldSet argument must either,
  /// - be empty, in which case it will be allocated and zeroed within applyAD, OR
  /// - be correctly allocated, i.e., should contain one Field for each Field in target, with each
  ///   source Field allocated to the size of the source FunctionSpace and the levels of the
  ///   corresponding target Field. In this case, target data is accumulated onto source.
  void applyAD(atlas::FieldSet & source, const atlas::FieldSet & target) const;

 private:
  void print(std::ostream &) const;

  const eckit::mpi::Comm & comm_;
  const atlas::FunctionSpace & source_fs_;
  const atlas::FunctionSpace & target_fs_;

  std::vector<std::vector<size_t>> mytarget_index_by_task_;
  std::vector<std::unique_ptr<LocalInterpolatorBase>> interp_;
  std::vector<int> mytarget_counts_;
  std::vector<int> mylocal_counts_;
};


// -----------------------------------------------------------------------------

}  // namespace oops
