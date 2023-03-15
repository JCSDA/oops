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

// TODO(Algo team): remove when eckit JEDI default version includes commit 2bda26c
#include <mpi.h>

#include <string>
#include <utility>
#include <vector>

#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"
// TODO(Algo team): remove when eckit JEDI default version includes commit 2bda26c
#include "eckit/mpi/Parallel.h"

#include "oops/util/Timer.h"

namespace util {
class DateTime;
}  // namespace util

namespace oops {
namespace mpi {

// ------------------------------------------------------------------------------------------------

/// Communicator with all MPI tasks (ie MPI_COMM_WORLD)
const eckit::mpi::Comm & world();

/// Communicator with each MPI task by itself
const eckit::mpi::Comm & myself();

// ------------------------------------------------------------------------------------------------

// TODO(Algo team): remove when eckit JEDI default version includes commit 2bda26c
MPI_Comm MPIComm(const eckit::mpi::Comm &comm);

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

template <typename SERIALIZABLE>
void sendReceiveReplace(const eckit::mpi::Comm & comm, SERIALIZABLE & sendrecvobj,
                        const int dest, const int sendtag, const int source, const int recvtag) {
  util::Timer timer("oops::mpi", "sendReceiveReplace");
  size_t sz = sendrecvobj.serialSize();
  std::vector<double> sendrecvbuf;
  sendrecvobj.serialize(sendrecvbuf);
  // TODO(Algo team): remove when eckit JEDI default version includes commit 2bda26c
  MPI_Comm mpicomm = MPIComm(comm);
  MPI_Status status;
  MPI_Sendrecv_replace(sendrecvbuf.data(), static_cast<int>(sz), MPI_DOUBLE,
                       dest, sendtag, source, recvtag, mpicomm, &status);
  // TODO(Algo team): uncomment when eckit JEDI default version includes commit 2bda26c
  // eckit::mpi::Status status = comm.sendReceiveReplace(sendrecvbuf.data(), sz,
  //                                                  dest, sendtag, source, recvtag);
  size_t ii = 0;
  sendrecvobj.deserialize(sendrecvbuf, ii);
  ASSERT(ii == sz);
}

// ------------------------------------------------------------------------------------------------

template <typename SERIALIZABLE>
void broadcast(const eckit::mpi::Comm & comm, SERIALIZABLE & obj, const size_t root) {
  util::Timer timer("oops::mpi", "broadcast");
  size_t sz = obj.serialSize();
  std::vector<double> buf;
  if (comm.rank() == root) {
    obj.serialize(buf);
  } else {
    buf.resize(sz);
  }
  comm.broadcast(buf, root);
  size_t ii = 0;
  obj.deserialize(buf, ii);
  ASSERT(ii == sz);
}

// ------------------------------------------------------------------------------------------------

template <typename SERIALIZABLE>
void reduceInPlace(const eckit::mpi::Comm & comm, SERIALIZABLE & obj, const size_t root) {
  util::Timer timer("oops::mpi", "reduceInPlace");
  size_t sz = obj.serialSize();
  std::vector<double> buf;
  obj.serialize(buf);
  // TODO(Algo team): remove when eckit JEDI default version includes commit 2bda26c
  MPI_Comm mpicomm = MPIComm(comm);
  if (comm.rank() == root) {
    MPI_Reduce(MPI_IN_PLACE, buf.data(), static_cast<int>(sz), MPI_DOUBLE, MPI_SUM,
               static_cast<int>(root), mpicomm);
  } else {
    MPI_Reduce(buf.data(), buf.data(), static_cast<int>(sz), MPI_DOUBLE, MPI_SUM,
               static_cast<int>(root), mpicomm);
  }
  // TODO(Algo team): uncomment when eckit JEDI default version includes commit 2bda26c
  // comm.reduceInPlace(buf.data(), sz, eckit::mpi::sum(), root);
  size_t ii = 0;
  if (comm.rank() == root) {
    obj.deserialize(buf, ii);
    ASSERT(ii == sz);
  }
}

// ------------------------------------------------------------------------------------------------

template <typename SERIALIZABLE>
void allReduceInPlace(const eckit::mpi::Comm & comm, SERIALIZABLE & obj) {
  util::Timer timer("oops::mpi", "allReduceInPlace");
  size_t sz = obj.serialSize();
  std::vector<double> buf;
  obj.serialize(buf);
  comm.allReduceInPlace(buf.data(), sz, eckit::mpi::sum());
  size_t ii = 0;
  obj.deserialize(buf, ii);
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

template <typename T>
void allGatherv(const eckit::mpi::Comm & comm, const std::vector<T> & send, std::vector<T> & recv) {
    eckit::mpi::Buffer<T> buffer(comm.size());
    comm.allGatherv(send.begin(), send.end(), buffer);
    recv = std::move(buffer.buffer);
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

/// \brief Perform the exclusive scan operation.
///
/// On output, `x` is set to the sum of the values of `x` passed to this function
/// on all ranks lower than the calling rank (and to 0 on rank 0).
void exclusiveScan(const eckit::mpi::Comm &comm, size_t &x);

}  // namespace mpi
}  // namespace oops
