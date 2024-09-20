/*
 * (C) Copyright 2023- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/generic/GlobalInterpolator.h"

#include <algorithm>
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

namespace {
  template<typename T>
  std::vector<T> flatten(const std::vector<std::vector<T>>& v) {
    size_t size = 0;
    for (const auto& inner : v) {
      size += inner.size();
    }
    std::vector<T> flat(size);
    size_t start = 0;
    for (auto const & inner : v) {
      std::copy(inner.begin(), inner.end(), flat.begin() + start);
      start += inner.size();
    }
    return flat;
  }

  template<typename T, typename IndexTy>
  std::vector<std::vector<T>> expand(const std::vector<T>& v, const std::vector<IndexTy>& sizes) {
    std::vector<std::vector<T>> expanded;
    expanded.reserve(sizes.size());
    IndexTy start = 0;
    for (const auto size : sizes) {
      expanded.emplace_back(v.begin() + start, v.begin() + start + size);
      start += size;
    }
    return expanded;
  }

  template<typename T>
  std::vector<T> scale(const std::vector<T>& v, T scale) {
    std::vector<T> scaled(v.size());
    std::transform(v.begin(), v.end(), scaled.begin(),
      [scale](const auto& x) { return x * scale; });
    return scaled;
  }

  std::vector<int> displsFromCounts(const std::vector<int>& counts) {
    std::vector<int> displs(counts.size());
    displs[0] = 0;
    std::partial_sum(counts.begin(), counts.end() - 1, displs.begin() + 1);
    return displs;
  }

  template<typename T>
  bool hasGivenSizes(const std::vector<std::vector<T>>& v, std::vector<int> test) {
    ASSERT(v.size() == test.size());
    auto testit = test.cbegin();
    return std::all_of(v.begin(), v.end(), [&testit](const auto& inner) {
        const auto is_same_size = inner.size() == static_cast<size_t>(*testit);
        ++testit;
        return is_same_size;
      });
  }

  template<typename T>
  std::vector<std::vector<T>> redistribute(
      const eckit::mpi::Comm& comm,
      const std::vector<std::vector<T>>& senddata,
      const std::vector<int>& sendcounts,
      const std::vector<int>& recvcounts,
      const int vec)
  {
    ASSERT(comm.size() == senddata.size() &&
           senddata.size() == sendcounts.size() &&
           sendcounts.size() == recvcounts.size());

    const auto sendcounts_v = scale(sendcounts, static_cast<int>(vec));
    ASSERT(hasGivenSizes(senddata, sendcounts_v));
    const auto recvcounts_v = scale(recvcounts, static_cast<int>(vec));
    const auto senddispls_v = displsFromCounts(sendcounts_v);
    const auto recvdispls_v = displsFromCounts(recvcounts_v);
    const auto recvcnt_v = recvdispls_v.back() + recvcounts_v.back();

    const auto sendbuf = flatten(senddata);
    std::vector<double> recvbuf(recvcnt_v);

    comm.allToAllv(sendbuf.data(), sendcounts_v.data(), senddispls_v.data(),
                   recvbuf.data(), recvcounts_v.data(), recvdispls_v.data());

    const auto recvdata = expand(recvbuf, recvcounts_v);
    return recvdata;
  }
}  // namespace


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

  // Infer mytarget/mylocal count info for use in alltoallv calls
  mytarget_counts_.resize(ntasks);
  mylocal_counts_.resize(ntasks);
  for (size_t i = 0; i < ntasks; ++i) {
    mytarget_counts_[i] = mytarget_latlon_by_task[i].size() / 2;
    mylocal_counts_[i] = mylocs_latlon_by_task[i].size() / 2;
  }

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
    nvars += (field.rank() == 1) ? 1 : field.shape(1);
    vars.push_back(field.name());
  }

  const size_t ntasks = comm_.size();

  source.haloExchange();

  std::vector<std::vector<double>> mylocal_interp(ntasks);
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    interp_[jtask]->apply(vars, source, mylocal_interp[jtask]);
  }

  // Gather results across MPI ranks
  const auto mytarget_interp = redistribute(comm_, mylocal_interp, mylocal_counts_,
                                            mytarget_counts_, nvars);

  // Copy data from vector<double> to atlas::FieldSet
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    const size_t nvals = mytarget_index_by_task_[jtask].size() * nvars;
    ASSERT(mytarget_interp[jtask].size() == nvals);
    if (nvals > 0) {
      LocalInterpolatorBase::bufferToFieldSet(vars, mytarget_index_by_task_[jtask],
                                              mytarget_interp[jtask], target);
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
    nvars += (field.rank() == 1) ? 1 : field.shape(1);
    vars.push_back(field.name());
  }

  const size_t ntasks = comm_.size();

  // (Adjoint of) Copy data from vector<double> to atlas::FieldSet
  std::vector<std::vector<double>> mytarget_interp(ntasks);
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    const size_t nvals = mytarget_index_by_task_[jtask].size() * nvars;
    mytarget_interp[jtask].resize(nvals, 0.0);
    if (nvals > 0) {
      LocalInterpolatorBase::bufferToFieldSetAD(vars, mytarget_index_by_task_[jtask],
                                                mytarget_interp[jtask], target);
    }
  }

  // (Adjoint of) Gather results across MPI ranks
  const auto mylocal_interp = redistribute(comm_, mytarget_interp, mytarget_counts_,
                                           mylocal_counts_, nvars);

  // (Adjoint of) Interpolate
  for (size_t jtask = 0; jtask < ntasks; ++jtask) {
    interp_[jtask]->applyAD(vars, source, mylocal_interp[jtask]);
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
