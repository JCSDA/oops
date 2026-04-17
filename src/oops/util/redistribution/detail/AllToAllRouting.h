/*
 * (C) Crown Copyright 2025-2026, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <vector>

#include "eckit/mpi/Comm.h"

namespace util {
namespace detail {

// ------------------------------------------------------------------------------------------------

/// \brief Helper class to facilitate allToAllv MPI operations.
class AllToAllRouting {
 public:
  AllToAllRouting() {}
  AllToAllRouting(const std::vector<int> sendCounts,
                  const std::vector<int> recvCounts);

  /// \brief For scaling the routing by the data size.
  AllToAllRouting& operator*=(const int rhs);

  /// \brief Execute the MPI allToAllv with given data buffers.
  template<typename T>
  void execute(const eckit::mpi::Comm& comm,
               const std::vector<T>& sendBuff,
               std::vector<T>& recvBuff) const {
    comm.allToAllv(sendBuff.data(), sendCounts_.data(), sendDispls_.data(),
                   recvBuff.data(), recvCounts_.data(), recvDispls_.data());
  }

  /// \brief Create a new routing with inverted receive and send arrays.
  AllToAllRouting invert() const;

 private:
  std::vector<int> sendCounts_;
  std::vector<int> sendDispls_;
  std::vector<int> recvCounts_;
  std::vector<int> recvDispls_;
};

/// \brief For scaling the routing by the data size.
AllToAllRouting operator*(const AllToAllRouting& lhs, const int rhs);

// ------------------------------------------------------------------------------------------------

}  // namespace detail
}  // namespace util
