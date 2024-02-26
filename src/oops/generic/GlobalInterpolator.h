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

  void apply(const atlas::FieldSet &, atlas::FieldSet &) const;
  void applyAD(atlas::FieldSet &, const atlas::FieldSet &) const;

 private:
  void print(std::ostream &) const;

  const eckit::mpi::Comm & comm_;
  std::vector<std::vector<size_t>> mytarget_index_by_task_;
  std::vector<std::unique_ptr<LocalInterpolatorBase>> interp_;
};


// -----------------------------------------------------------------------------

}  // namespace oops
