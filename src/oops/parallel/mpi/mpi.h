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

    const eckit::mpi::Comm& comm();

//-------------------------------------------------------------------------------------------------
// allGather for eigen vectors
//-------------------------------------------------------------------------------------------------
    void allGather(const Eigen::VectorXd &, std::vector<Eigen::VectorXd> &);

//-------------------------------------------------------------------------------------------------
// Send and receive objects with a serialize/deserialize method
//-------------------------------------------------------------------------------------------------
    // Non blocking DO NOT USE
    template<typename TYPE>
    eckit::mpi::Request iSend(const TYPE & sendbuf, int & source, int & tag) {
      std::vector<double> v_sendbuf;
      sendbuf.serialize(v_sendbuf);
      return(comm().iSend(v_sendbuf.data(), v_sendbuf.size(), source, tag));
//    v_sendbuf gets deallocated here, potentially before it has been sent
      ABORT("Bug in iSend");
    }

    template<typename TYPE>
    eckit::mpi::Request iReceive(TYPE & to_fill, int & source, int & tag) {
      size_t sz = to_fill.serialSize();
      std::vector<double> v_recv(sz);
      eckit::mpi::Request request = comm().iReceive(v_recv.data(), sz, source, tag);
      ABORT("Bug in iReceive");
//    v_recv has not been received at this point, it will be after wait
      to_fill.deserialize(v_recv);
      return request;
    }

    // Blocking
    template<typename TYPE>
    void send(const TYPE & sendbuf, int & source, int & tag) {
      std::vector<double> v_sendbuf;
      sendbuf.serialize(v_sendbuf);
      comm().send(v_sendbuf.data(), v_sendbuf.size(), source, tag);
    }

    template<typename TYPE>
    void receive(TYPE & to_fill, int & source, int & tag) {
      size_t sz = to_fill.serialSize();
      std::vector<double> v_recv(sz);
      eckit::mpi::Status status = comm().receive(v_recv.data(), sz, source, tag);
      to_fill.deserialize(v_recv);
    }

  }  // namespace mpi
}  // namespace oops
