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

#include <vector>

#include "eckit/mpi/Comm.h"
#include "oops/util/abor1_cpp.h"

namespace oops {
  namespace mpi {

/// Default communicator with all MPI tasks (ie MPI_COMM_WORLD)
    const eckit::mpi::Comm & world();

/// Default communicator with each MPI task by itself
    const eckit::mpi::Comm & myself();

//-------------------------------------------------------------------------------------------------
// allGather for eigen vectors
//-------------------------------------------------------------------------------------------------
    void allGather(const Eigen::VectorXd &, std::vector<Eigen::VectorXd> &);

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

//-------------------------------------------------------------------------------------------------
// Send and receive objects with a serialize/deserialize method
//-------------------------------------------------------------------------------------------------
    // Non blocking DO NOT USE
    template<typename TYPE>
    eckit::mpi::Request iSend(const eckit::mpi::Comm & comm,
                              const TYPE & sendbuf, int & source, int & tag) {
      std::vector<double> v_sendbuf;
      sendbuf.serialize(v_sendbuf);
      return(comm.iSend(v_sendbuf.data(), v_sendbuf.size(), source, tag));
//    v_sendbuf gets deallocated here, potentially before it has been sent
      ABORT("Bug in iSend");
    }

    template<typename TYPE>
    eckit::mpi::Request iReceive(const eckit::mpi::Comm & comm,
                                 TYPE & to_fill, int & source, int & tag) {
      size_t sz = to_fill.serialSize();
      std::vector<double> v_recv(sz);
      eckit::mpi::Request request = comm.iReceive(v_recv.data(), sz, source, tag);
      ABORT("Bug in iReceive");
//    v_recv has not been received at this point, it will be after wait
      to_fill.deserialize(v_recv);
      return request;
    }

    // Blocking
    template<typename TYPE>
    void send(const eckit::mpi::Comm & comm, const TYPE & sendbuf, int & source, int & tag) {
      std::vector<double> v_sendbuf;
      sendbuf.serialize(v_sendbuf);
      comm.send(v_sendbuf.data(), v_sendbuf.size(), source, tag);
    }

    template<typename TYPE>
    void receive(const eckit::mpi::Comm & comm, TYPE & to_fill, int & source, int & tag) {
      size_t sz = to_fill.serialSize();
      std::vector<double> v_recv(sz);
      eckit::mpi::Status status = comm.receive(v_recv.data(), sz, source, tag);
      to_fill.deserialize(v_recv);
    }

  }  // namespace mpi
}  // namespace oops
