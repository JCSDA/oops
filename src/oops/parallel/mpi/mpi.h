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

#include <vector>

#include "eckit/mpi/Comm.h"

namespace oops {
  namespace mpi {

    const eckit::mpi::Comm& comm();

//-------------------------------------------------------------------------------------------------
// Send a receive any object with serialize/deserialize methods
//-------------------------------------------------------------------------------------------------
    // Non blocking
    template<typename TYPE>
    eckit::mpi::Request iSend(const TYPE & sendbuf, int & source, int & tag) {
      std::vector<double> v_sendbuf = sendbuf.serialize();
      return(comm().iSend(v_sendbuf.data(), v_sendbuf.size(), source, tag));
    }

    template<typename TYPE>
    eckit::mpi::Request iReceive(std::vector<double> & v_recv, TYPE & to_fill,
                                 int & source, int & tag) {
      eckit::mpi::Request request = comm().iReceive(v_recv.data(), v_recv.size(), source, tag);
      to_fill.deserialize(v_recv);
      return request;
    }

    // Blocking
    template<typename TYPE>
    void send(const TYPE & sendbuf, int & source, int & tag) {
      std::vector<double> v_sendbuf = sendbuf.serialize();
      comm().send(v_sendbuf.data(), v_sendbuf.size(), source, tag);
    }

    template<typename TYPE>
    void receive(std::vector<double> & v_recv, TYPE & to_fill, int & source,
                 int & tag) {
      eckit::mpi::Status status = comm().receive(v_recv.data(), v_recv.size(), source, tag);
      to_fill.deserialize(v_recv);
    }

  }  // namespace mpi
}  // namespace oops
