// (C) Crown Copyright 2026 Met Office
//
// This software is licensed under the terms of the Apache Licence Version 2.0
// which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.

#pragma once

#include <cstddef>
#include <utility>
#include <vector>

#include "eckit/exception/Exceptions.h"

namespace oops {

/// \brief Describes how observation data for one MPI task is laid out in a flat buffer.
///
/// GetValues packs interpolation results into one flat buffer per MPI task. Within a
/// task's buffer the data is arranged as one contiguous block per sampling method, in
/// ascending sampling-method order; within a block the ordering is
///   variable (slowest) -> observation (middle) -> level (fastest).
/// A block therefore holds `obsCount(jtask, jsm) * varsize(jsm)` values, where
/// `varsize(jsm)` is the sum over that sampling method's variables of their number of
/// levels.
///
/// The same layout class is used for two different groupings of observations:
///  - grouped by obs-owning task ("ot"), for the buffers that are sent, and
///  - grouped by model-owning / interpolating task ("mt"), for the buffers that are
///    received.
///
/// The number of levels per variable is not stored here, because a single grouping of
/// observations is used with two different variable sets (the GeoVaLs variables on the
/// nonlinear path, the linear variables on the TL/AD path). It is passed in per query.
class ObsBufferDescriptor {
 public:
  ObsBufferDescriptor() = default;

  /// \param obsCounts obsCounts[jtask][jsm] = number of observations held by task
  ///                  `jtask` for sampling method `jsm`. Must be rectangular.
  explicit ObsBufferDescriptor(std::vector<std::vector<size_t>> obsCounts)
    : obsCounts_(std::move(obsCounts)) {
    for (const auto & perTask : obsCounts_) {
      ASSERT(perTask.size() == nsms());
    }
  }

  size_t ntasks() const { return obsCounts_.size(); }
  size_t nsms() const { return obsCounts_.empty() ? 0 : obsCounts_[0].size(); }

  size_t obsCount(const size_t jtask, const size_t jsm) const {
    return obsCounts_[jtask][jsm];
  }

  /// \brief Total number of values in task `jtask`'s buffer, across all sampling methods.
  /// \param varsizes_sm varsizes_sm[jsm] = total number of levels summed over the
  ///                    variables of sampling method `jsm`.
  size_t bufferSize(const size_t jtask, const std::vector<size_t> & varsizes_sm) const {
    ASSERT(varsizes_sm.size() == nsms());
    size_t size = 0;
    for (size_t jsm = 0; jsm < nsms(); ++jsm) {
      size += obsCounts_[jtask][jsm] * varsizes_sm[jsm];
    }
    return size;
  }

 private:
  /// obsCounts_[jtask][jsm]
  std::vector<std::vector<size_t>> obsCounts_;
};

}  // namespace oops
