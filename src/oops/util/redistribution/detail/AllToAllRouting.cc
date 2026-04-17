/*
 * (C) Crown Copyright 2025-2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "eckit/exception/Exceptions.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/redistribution/detail/AllToAllRouting.h"

namespace util {
namespace detail {

// ------------------------------------------------------------------------------------------------

AllToAllRouting::AllToAllRouting(const std::vector<int> sendCounts,
                                 const std::vector<int> recvCounts) :
  sendCounts_(sendCounts), recvCounts_(recvCounts),
  sendDispls_(sendCounts.size(), 0),
  recvDispls_(recvCounts.size(), 0)
{
  ASSERT(sendCounts_.size() == recvCounts_.size());

  // Calculate displacement arrays.
  // E.g, if sendCounts = {1, 2, 3, 4};
  //         sendDispls = {0, 1, 3, 6};    (cumulative sum)
  const auto calcDisp = [](const std::vector<int>& counts) {
    std::vector<int> disp;
    disp.reserve(counts.size());

    size_t rollingTotal = 0;

    disp.emplace_back(rollingTotal);
    for (size_t idx = 0; idx < counts.size() - 1; ++idx) {
      rollingTotal += counts[idx];
      disp.emplace_back(rollingTotal);
    }
    return disp;
  };
  sendDispls_ = calcDisp(sendCounts_);
  recvDispls_ = calcDisp(recvCounts_);
}

// ------------------------------------------------------------------------------------------------

AllToAllRouting& AllToAllRouting::operator*=(const int rhs) {
  for (size_t i = 0; i < sendCounts_.size(); ++i) {
    sendCounts_[i] *= rhs;
    sendDispls_[i] *= rhs;
    recvCounts_[i] *= rhs;
    recvDispls_[i] *= rhs;
  }
  return *this;
}

// ------------------------------------------------------------------------------------------------

AllToAllRouting AllToAllRouting::invert() const {
  return AllToAllRouting(recvCounts_, sendCounts_);
}

// ------------------------------------------------------------------------------------------------

AllToAllRouting operator*(const AllToAllRouting& lhs, const int rhs) {
  AllToAllRouting result(lhs);
  return result *= rhs;
}

// ------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace util

