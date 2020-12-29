/*
 * (C) Copyright 2013 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "oops/mpi/mpi.h"

#include <string>

#include "eckit/exception/Exceptions.h"

namespace oops {
namespace mpi {

// ------------------------------------------------------------------------------------------------

const eckit::mpi::Comm & world() {
  return eckit::mpi::comm();
}

// ------------------------------------------------------------------------------------------------

const eckit::mpi::Comm & myself() {
  return eckit::mpi::self();
}

// ------------------------------------------------------------------------------------------------

void gather(const eckit::mpi::Comm & comm, const std::vector<double> & send,
            std::vector<double> & recv, const size_t root) {
  size_t ntasks = comm.size();
  if (ntasks > 1) {
    int mysize = send.size();
    std::vector<int> sizes(ntasks);
    comm.allGather(mysize, sizes.begin(), sizes.end());
    std::vector<int> displs(ntasks);
    size_t rcvsz = sizes[0];
    displs[0] = 0;
    for (size_t jj = 1; jj < ntasks; ++jj) {
      displs[jj] = displs[jj - 1] + sizes[jj - 1];
      rcvsz += sizes[jj];
    }
    if (comm.rank() == root) recv.resize(rcvsz);

    comm.gatherv(send, recv, sizes, displs, root);
  } else {
    recv = send;
  }
}

// ------------------------------------------------------------------------------------------------

void allGather(const eckit::mpi::Comm & comm,
               const Eigen::VectorXd & sendbuf, Eigen::MatrixXd & recvbuf) {
  const int ntasks = comm.size();
  int buf_size = sendbuf.size();

  std::vector<double> vbuf(sendbuf.data(), sendbuf.data() + buf_size);
  std::vector<double> vbuf_total(ntasks * buf_size);

  std::vector<int> recvcounts(ntasks);
  for (int ii = 0; ii < ntasks; ++ii) recvcounts[ii] = buf_size;

  std::vector<int> displs(ntasks);
  for (int ii = 0; ii < ntasks; ++ii) displs[ii] = ii * buf_size;

  comm.allGatherv(vbuf.begin(), vbuf.end(),
                  vbuf_total.begin(), recvcounts.data(), displs.data());

  for (int ii = 0; ii < ntasks; ++ii) {
    std::vector<double> vloc(vbuf_total.begin() + ii * buf_size,
                             vbuf_total.begin() + (ii + 1) * buf_size);
    Eigen::VectorXd my_vect = Eigen::Map<Eigen::VectorXd, Eigen::Unaligned>(vloc.data(),
                                                                            vloc.size());
    recvbuf.col(ii) = my_vect;
  }
}

// ------------------------------------------------------------------------------------------------

}  // namespace mpi
}  // namespace oops
