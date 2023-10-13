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

/// Communicator with all MPI tasks (ie MPI_COMM_WORLD)
const eckit::mpi::Comm & world();

/// Communicator with each MPI task by itself
const eckit::mpi::Comm & myself();

const eckit::mpi::Comm & clone(const eckit::mpi::Comm &);

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
  eckit::mpi::Status status = comm.sendReceiveReplace(sendrecvbuf.data(), sz,
                                                      dest, sendtag, source, recvtag);
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
  comm.reduceInPlace(buf.data(), sz, eckit::mpi::sum(), root);
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

// ------------------------------------------------------------------------------------------------
// MPI broadcast utilities based on eckit broadcast.

/// \brief broadcast a vector variable via the eckit broadcast
/// \param comm eckit communicator group
/// \param vectorVar vector for broadcasting
/// \param root root rank for broadcasting
template <typename VecType>
void broadcastVector(const eckit::mpi::Comm & comm, std::vector<VecType> & vectorVar,
                     const size_t root) {
    // eckit broadcast support vectors, but you need to have the vectors identically
    // sized on both sides before doing the broadcast. This routine will broadcast
    // the vector size so the receiving end can resize properly.
    int vecSize;
    if (comm.rank() == root) {
        vecSize = vectorVar.size();
        comm.broadcast(vecSize, root);
        comm.broadcast(vectorVar, root);
    } else {
        comm.broadcast(vecSize, root);
        vectorVar.resize(vecSize);
        comm.broadcast(vectorVar, root);
    }
}

/// @brief broadcast a bool variable via the eckit broadcast
/// @param comm eckit communicator group
/// @param boolVar variable for broadcasting
/// @param root root rank for broadcasting
void broadcastBool(const eckit::mpi::Comm & comm, bool & boolVar, const size_t root);

/// \brief broadcast a string variable via the eckit broadcast
/// \param comm eckit communicator group
/// \param stringVar string for broadcasting
/// \param root root rank for broadcasting
void broadcastString(const eckit::mpi::Comm & comm, std::string & stringVar, const size_t root);

// ------------------------------------------------------------------------------------------------
// MPI send/receive utilities based on eckit send/receive.

/// \brief send a string variable via the eckit send
/// \param comm eckit communicator group
/// \param stringVar string for sending
/// \param toRank rank we are sending to
void sendString(const eckit::mpi::Comm & comm, const std::string & stringVar, const int toRank);

/// \brief receive a string variable via the eckit receive
/// \param comm eckit communicator group
/// \param stringVar string for receiving
/// \param fromRank rank we are receiving from
void receiveString(const eckit::mpi::Comm & comm, std::string & stringVar, const int fromRank);

}  // namespace mpi
}  // namespace oops
