/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#pragma once

#include <Eigen/Dense>

#include <string>
#include <utility>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/Timer.h"

namespace util {
class DateTime;
}  // namespace util

namespace oops {
namespace mpi {

// ------------------------------------------------------------------------------------------------

/// Default communicator with all MPI tasks (ie MPI_COMM_WORLD)
const eckit::mpi::Comm & world();

/// Default communicator with each MPI task by itself
const eckit::mpi::Comm & myself();

// ------------------------------------------------------------------------------------------------

/// Extend eckit Comm for Serializable oops objects

template <typename SERIALIZABLE>
void send(const eckit::mpi::Comm & comm, const SERIALIZABLE & sendobj,
          const int dest, const int tag) {
  util::Timer timer("oops::mpi", "send");
  std::vector<double> sendbuf;
  sendobj.serialize(sendbuf);
  comm.send(sendbuf.data(), sendbuf.size(), dest, tag);
}

// ------------------------------------------------------------------------------------------------

template <typename SERIALIZABLE>
void receive(const eckit::mpi::Comm & comm, SERIALIZABLE & recvobj,
             const int source, const int tag) {
  util::Timer timer("oops::mpi", "receive");
  size_t sz = recvobj.serialSize();
  std::vector<double> recvbuf(sz);
  eckit::mpi::Status status = comm.receive(recvbuf.data(), sz, source, tag);
  size_t ii = 0;
  recvobj.deserialize(recvbuf, ii);
  ASSERT(ii == sz);
}

// ------------------------------------------------------------------------------------------------

void gather(const eckit::mpi::Comm & comm, const std::vector<double> & send,
            std::vector<double> & recv, const size_t root);

// ------------------------------------------------------------------------------------------------

template <typename SERIALIZABLE>
void gather(const eckit::mpi::Comm & comm, const std::vector<SERIALIZABLE> & send,
            std::vector<SERIALIZABLE> & recv, const size_t root) {
  if (comm.size() > 1) {
    std::vector<double> sendbuf;
    std::vector<double> recvbuf;

    for (const SERIALIZABLE & jsend : send) jsend.serialize(sendbuf);

    gather(comm, sendbuf, recvbuf, root);

    if (comm.rank() == root) {
      size_t indx = 0;
      for (SERIALIZABLE & jrecv : recv) jrecv.deserialize(recvbuf, indx);
    }
  } else {
    recv = send;
  }
}

// ------------------------------------------------------------------------------------------------
// allGather for eigen vectors
// ------------------------------------------------------------------------------------------------
void allGather(const eckit::mpi::Comm & comm,
               const Eigen::VectorXd &, Eigen::MatrixXd &);

// ------------------------------------------------------------------------------------------------

/// \brief A wrapper around the MPI *all gather* operation for serializable types.
///
/// The *all gather* operation gathers data from all tasks and delivers the combined data to all
/// tasks. This wrapper performs that operation for collections of non-primitive types that
/// nevertheless support the OOPS serialization interface, i.e. provide the functions
///
///     void serialize(std::vector<double> &vect) const;
///     void deserialize(const std::vector<double> &vect, size_t &current);
///
/// An example of such a type is util::DateTime.
///
/// \param comm
///   Communicator.
/// \param first, last
///   Range of values to be delivered from this task to all other tasks.
/// \param recvbuf
///   Output iterator to the beginning of the range to receive the combined data from all tasks.
template <typename CIter, typename Iter>
void allGathervUsingSerialize(const eckit::mpi::Comm &comm, CIter first, CIter last,
                              Iter recvbuf) {
  std::vector<double> serializedLocalData;
  for (CIter it = first; it != last; ++it)
    it->serialize(serializedLocalData);

  eckit::mpi::Buffer<double> buffer(comm.size());
  comm.allGatherv(serializedLocalData.begin(), serializedLocalData.end(), buffer);

  size_t numDeserializedDoubles = 0;
  for (Iter it = recvbuf; numDeserializedDoubles != buffer.buffer.size(); ++it)
    it->deserialize(buffer.buffer, numDeserializedDoubles);
}

// ------------------------------------------------------------------------------------------------

// The following functions simplify the allGatherv operation on vectors (reducing it to a single
// function call).

// ------------------------------------------------------------------------------------------------

/// \brief Perform the MPI *all gather* operation on a vector of "plain old data".
///
/// This operation gathers data from all tasks and delivers the combined data to all tasks.
///
/// \tparam T must be a type for which there exists a specialization of eckit::mpi::Data::Type.
///
/// \param[in] comm
///   Communicator.
/// \param[inout] x
///   On input, data owned by this task that need to be delivered to all other tasks. On output,
///   combined data received from all tasks (concatenated in the order of increasing task ranks).
template <typename T>
void allGatherv(const eckit::mpi::Comm & comm, std::vector<T> &x) {
    eckit::mpi::Buffer<T> buffer(comm.size());
    comm.allGatherv(x.begin(), x.end(), buffer);
    x = std::move(buffer.buffer);
}

// ------------------------------------------------------------------------------------------------

/// \brief Perform the MPI *all gather* operation on a vector of DateTime objects.
///
/// This operation gathers data from all tasks and delivers the combined data to all tasks.
///
/// \param[in] comm
///   Communicator.
/// \param[inout] x
///   On input, data owned by this task that need to be delivered to all other tasks. On output,
///   combined data received from all tasks (concatenated in the order of increasing task ranks).
void allGatherv(const eckit::mpi::Comm & comm, std::vector<util::DateTime> &x);

// ------------------------------------------------------------------------------------------------

/// \brief Perform the MPI *all gather* operation on a vector of DateTime objects.
///
/// This operation gathers data from all tasks and delivers the combined data to all tasks.
///
/// \param[in] comm
///   Communicator.
/// \param[inout] x
///   On input, data owned by this task that need to be delivered to all other tasks. On output,
///   combined data received from all tasks (concatenated in the order of increasing task ranks).
void allGatherv(const eckit::mpi::Comm & comm, std::vector<std::string> &x);

// ------------------------------------------------------------------------------------------------

}  // namespace mpi
}  // namespace oops
