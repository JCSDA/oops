/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/util/Geometry.h"
#include "atlas/util/KDTree.h"

#include "eckit/config/Configuration.h"

#include "oops/base/Variables.h"
#include "oops/generic/UnstructuredInterpolator.h"
#include "oops/interface/LocalInterpolator.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"


namespace oops {

// -----------------------------------------------------------------------------

// Select UnstructuredInterpolator (default) or LocalInterpolator.
template<typename MODEL, typename = void>
struct SelectInterp {
  typedef UnstructuredInterpolator type;
};

template<typename MODEL>
struct SelectInterp<MODEL, cpp17::void_t<typename MODEL::LocalInterpolator>> {
  typedef LocalInterpolator<MODEL> type;
};


template<typename MODEL>
class GlobalInterpolator : public util::Printable {
  typedef Geometry<MODEL>  Geometry_;
  typedef typename SelectInterp<MODEL>::type Interp_;

 public:
  // Constructor used when interpolating to random locations for testing
  GlobalInterpolator(const eckit::Configuration &,
                     const Geometry_ &,
                     const atlas::FunctionSpace &,
                     const eckit::mpi::Comm &);
  ~GlobalInterpolator() = default;

  void apply(const atlas::FieldSet &, atlas::FieldSet &) const;
  void applyAD(atlas::FieldSet &, const atlas::FieldSet &) const;

 private:
  void setupInterpolators(const eckit::Configuration &,
                          const Geometry_ &,
                          const std::vector<double>&,
                          const std::vector<double>&);

  void print(std::ostream &) const;

  const eckit::mpi::Comm & comm_;
  std::vector<std::vector<size_t>> mytarget_index_by_task_;
  std::vector<std::unique_ptr<Interp_>> interp_;
};

// -----------------------------------------------------------------------------

template <typename MODEL>
GlobalInterpolator<MODEL>::GlobalInterpolator(
    const eckit::Configuration & config,
    const Geometry_ & source_grid,
    const atlas::FunctionSpace & target_fs,
    const eckit::mpi::Comm & comm)
  : comm_(comm)
{
  Log::trace() << "GlobalInterpolator::GlobalInterpolator start" << std::endl;

  // Extract target coords
  const auto lonlat = atlas::array::make_view<double, 2>(target_fs.lonlat());
  const size_t npts = target_fs->size();
  std::vector<double> target_lats(npts);
  std::vector<double> target_lons(npts);
  for (size_t jj = 0; jj < npts; ++jj) {
    target_lats[jj] = lonlat(jj, 1);
    target_lons[jj] = lonlat(jj, 0);
    if (target_lons[jj] < 0.0) target_lons[jj] += 360.0;
  }

  setupInterpolators(config, source_grid, target_lats, target_lons);

  Log::trace() << "GlobalInterpolator::GlobalInterpolator done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void GlobalInterpolator<MODEL>::setupInterpolators(
    const eckit::Configuration & config,
    const Geometry_ & source_grid,
    const std::vector<double> & target_lats,
    const std::vector<double> & target_lons)
{
  const size_t ntasks = comm_.size();
  const size_t npts = target_lats.size();

  // Exchange target coords, find targets to be interpolated on this task
  mytarget_index_by_task_.resize(ntasks);
  std::vector<std::vector<double>> mytarget_latlon_by_task(ntasks);
  for (size_t jt = 0; jt < npts; ++jt) {
    const size_t itask = source_grid.closestTask(target_lats[jt], target_lons[jt]);
    mytarget_index_by_task_[itask].push_back(jt);
    mytarget_latlon_by_task[itask].push_back(target_lats[jt]);
    mytarget_latlon_by_task[itask].push_back(target_lons[jt]);
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
    interp_[jtask].reset(new Interp_(config, source_grid, lats, lons));
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void GlobalInterpolator<MODEL>::apply(const atlas::FieldSet & source,
                                      atlas::FieldSet & target) const {
  ASSERT(source.field_names() == target.field_names());
  for (const auto & src_field : source) {
    const auto & tgt_field = target.field(src_field.name());
    ASSERT(src_field.rank() == tgt_field.rank());
    ASSERT(src_field.shape(1) == tgt_field.shape(1));
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
      Interp_::bufferToFieldSet(vars, mytarget_index_by_task_[jtask], recvinterp[jtask], target);
    }
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void GlobalInterpolator<MODEL>::applyAD(atlas::FieldSet & source,
                                        const atlas::FieldSet & target) const {
  ASSERT(target.field_names() == source.field_names());
  for (const auto & tgt_field : target) {
    const auto & src_field = source.field(tgt_field.name());
    ASSERT(tgt_field.rank() == src_field.rank());
    ASSERT(tgt_field.shape(1) == src_field.shape(1));
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
      Interp_::bufferToFieldSetAD(vars, mytarget_index_by_task_[jtask], recvinterp[jtask], target);
    }
  }

  // (Adjoint of) Gather results across MPI ranks
  std::vector<std::vector<double>> locinterp(ntasks);
  comm_.allToAll(recvinterp, locinterp);

  // (Adjoint of) Interpolate
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    interp_[jtask]->applyAD(vars, source, locinterp[jtask]);
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void GlobalInterpolator<MODEL>::print(std::ostream & os) const
{
  os << "GlobalInterpolator<" << MODEL::name() << ">";
}

// -----------------------------------------------------------------------------

}  // namespace oops
