/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/generic/GlobalInterpolator.h"

#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/util/Geometry.h"

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/GeometryData.h"
#include "oops/base/Variables.h"
#include "oops/generic/AtlasInterpolator.h"
#include "oops/generic/UnstructuredInterpolator.h"
#include "oops/util/Logger.h"


namespace oops {

// -----------------------------------------------------------------------------

GlobalInterpolator::GlobalInterpolator(
    const eckit::Configuration & config,
    const GeometryData & source_grid,
    const atlas::FunctionSpace & target_fs,
    const eckit::mpi::Comm & comm)
  : comm_(comm), source_fs_(source_grid.functionSpace()), target_fs_(target_fs)
{
  Log::trace() << "GlobalInterpolator::GlobalInterpolator start" << std::endl;

  const std::string interp_type = config.getString("local interpolator type");
  if (interp_type != "atlas interpolator"
      && interp_type != "oops unstructured grid interpolator") {
    throw eckit::BadParameter("Unknown local interpolator type: " + interp_type);
  }

  const size_t ntasks = comm_.size();
  mytarget_index_by_task_.resize(ntasks);
  std::vector<std::vector<double>> mytarget_latlon_by_task(ntasks);

  // Extract target coordinates, find which task will be responsible for interpolating to each one
  const auto lonlat = atlas::array::make_view<double, 2>(target_fs_.lonlat());
  const auto ghost = atlas::array::make_view<int, 1>(target_fs_.ghost());
  for (size_t jj = 0; jj < lonlat.shape(0); ++jj) {
    if (ghost(jj) == 0) {
      double lon = lonlat(jj, 0);
      if (lon < 0.0) lon += 360.0;
      const double lat = lonlat(jj, 1);
      const size_t itask = source_grid.closestTask(lat, lon);
      mytarget_index_by_task_[itask].push_back(jj);
      mytarget_latlon_by_task[itask].push_back(lat);
      mytarget_latlon_by_task[itask].push_back(lon);
    }
  }

  std::vector<std::vector<double>> mylocs_latlon_by_task(ntasks);
  comm_.allToAll(mytarget_latlon_by_task, mylocs_latlon_by_task);

  interp_.resize(ntasks);
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    const size_t ntargets = mylocs_latlon_by_task[jtask].size() / 2;
    std::vector<double> lats(ntargets);
    std::vector<double> lons(ntargets);
    size_t ii = 0;
    for (size_t jt = 0; jt < ntargets; ++jt) {
      lats[jt] = mylocs_latlon_by_task[jtask][ii];
      lons[jt] = mylocs_latlon_by_task[jtask][ii + 1];
      ii += 2;
    }
    ASSERT(mylocs_latlon_by_task[jtask].size() == ii);

    if (interp_type == "atlas interpolator") {
      interp_[jtask].reset(new AtlasInterpolator(config, source_grid, lats, lons));
    } else if (interp_type == "oops unstructured grid interpolator") {
      interp_[jtask].reset(new UnstructuredInterpolator(config, source_grid, lats, lons));
    }
  }

  Log::trace() << "GlobalInterpolator::GlobalInterpolator done" << std::endl;
}

// -----------------------------------------------------------------------------

void GlobalInterpolator::apply(const atlas::FieldSet & source,
                               atlas::FieldSet & target) const {
  if (target.empty()) {
    // Allocate target Fields to match source
    for (const auto & src_field : source) {
      atlas::Field tgt_field = target_fs_.createField<double>(
          atlas::option::name(src_field.name()) | atlas::option::levels(src_field.levels()));
      tgt_field.metadata() = src_field.metadata();
      target.add(tgt_field);
    }
  } else {
    // Verify target Fields correctly allocated
    ASSERT(target.field_names() == source.field_names());
    for (const auto & src_field : source) {
      const auto & tgt_field = target.field(src_field.name());
      ASSERT(tgt_field.rank() == src_field.rank());
      ASSERT(tgt_field.shape(0) == target_fs_.size());
      ASSERT(tgt_field.shape(1) == src_field.shape(1));
    }
  }

  // For a simple interface, construct variables from all fields in fieldset
  size_t nvars = 0;
  oops::Variables vars{};
  for (const auto & field : source) {
    const size_t rank = field.rank();
    if (rank == 1) {
      nvars += 1;
    } else {
      nvars += field.shape(1);
    }
    vars.push_back(field.name());
  }

  const size_t ntasks = comm_.size();

  source.haloExchange();

  std::vector<std::vector<double>> locinterp(ntasks);
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    interp_[jtask]->apply(vars, source, locinterp[jtask]);
  }

  // Gather results across MPI ranks
  std::vector<std::vector<double>> recvinterp(ntasks);
  comm_.allToAll(locinterp, recvinterp);

  // Copy data from vector<double> to atlas::FieldSet
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    const size_t nvals = mytarget_index_by_task_[jtask].size() * nvars;
    ASSERT(recvinterp[jtask].size() == nvals);
    if (nvals > 0) {
      LocalInterpolatorBase::bufferToFieldSet(vars, mytarget_index_by_task_[jtask],
                                              recvinterp[jtask], target);
    }
  }

  target.set_dirty();
}

// -----------------------------------------------------------------------------

void GlobalInterpolator::applyAD(atlas::FieldSet & source,
                                 const atlas::FieldSet & target) const {
  if (source.empty()) {
    // Allocate source Fields to match target, and zero before adjoint accumulation
    for (const auto & tgt_field : target) {
      atlas::Field src_field = source_fs_.createField<double>(
          atlas::option::name(tgt_field.name()) | atlas::option::levels(tgt_field.levels()));
      src_field.metadata() = tgt_field.metadata();
      auto src_view = atlas::array::make_view<double, 2>(src_field);
      src_view.assign(0.0);
      source.add(src_field);
    }
  } else {
    // Verify source Fields correctly allocated
    ASSERT(source.field_names() == target.field_names());
    for (const auto & tgt_field : target) {
      const auto & src_field = source.field(tgt_field.name());
      ASSERT(src_field.rank() == tgt_field.rank());
      ASSERT(src_field.shape(0) == source_fs_.size());
      ASSERT(src_field.shape(1) == tgt_field.shape(1));
    }
  }

  // For a simple interface, construct variables from all fields in fieldset
  size_t nvars = 0;
  oops::Variables vars{};
  for (const auto & field : target) {
    const size_t rank = field.rank();
    if (rank == 1) {
      nvars += 1;
    } else {
      nvars += field.shape(1);
    }
    vars.push_back(field.name());
  }

  const size_t ntasks = comm_.size();

  // (Adjoint of) Copy data from vector<double> to atlas::FieldSet
  std::vector<std::vector<double>> recvinterp(ntasks);
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    const size_t nvals = mytarget_index_by_task_[jtask].size() * nvars;
    recvinterp[jtask].resize(nvals, 0.0);
    if (nvals > 0) {
      LocalInterpolatorBase::bufferToFieldSetAD(vars, mytarget_index_by_task_[jtask],
                                                recvinterp[jtask], target);
    }
  }

  // (Adjoint of) Gather results across MPI ranks
  std::vector<std::vector<double>> locinterp(ntasks);
  comm_.allToAll(recvinterp, locinterp);

  // (Adjoint of) Interpolate
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    interp_[jtask]->applyAD(vars, source, locinterp[jtask]);
  }

  source.adjointHaloExchange();
  source.set_dirty();
}

// -----------------------------------------------------------------------------

void GlobalInterpolator::print(std::ostream & os) const
{
  os << "GlobalInterpolator";
}

// -----------------------------------------------------------------------------

}  // namespace oops
