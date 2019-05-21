/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "oops/parallel/mpi/mpi.h"

namespace oops {
  namespace mpi {

    const eckit::mpi::Comm& comm() {
      return eckit::mpi::comm();
    }

    void allGather(const Eigen::VectorXd & sendbuf, std::vector<Eigen::VectorXd> & recvbuf) {
      const int ntasks = comm().size();
      int buf_size = sendbuf.size();

      std::vector<double> vbuf(sendbuf.data(), sendbuf.data() + buf_size);
      std::vector<double> vbuf_total(ntasks * buf_size);

      std::vector<int> recvcounts(ntasks);
      for (int ii = 0; ii < ntasks; ++ii) recvcounts[ii] = buf_size;

      std::vector<int> displs(ntasks);
      for (int ii = 0; ii < ntasks; ++ii) displs[ii] = ii * buf_size;

      comm().allGatherv(vbuf.begin(), vbuf.end(),
                        vbuf_total.begin(), recvcounts.data(), displs.data());

      for (int ii = 0; ii < ntasks; ++ii) {
        std::vector<double> vloc(vbuf_total.begin() + ii * buf_size,
                                 vbuf_total.begin() + (ii + 1) * buf_size);
        Eigen::VectorXd my_vect = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(vloc.data(),
                                                                                vloc.size());
        recvbuf[ii] = my_vect;
      }
    }
  }  // namespace mpi
}  // namespace oops
