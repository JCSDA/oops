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

#include <numeric>  // for accumulate()
#include <string>
#include <utility>

#include "eckit/exception/Exceptions.h"
#include "oops/util/DateTime.h"

namespace {

// Helper functions used by the implementation of the specialization of allGatherv for a vectors
// of strings

/// \brief Join strings into a single character array before MPI transfer.
///
/// \param strings
///   Strings to join.
///
/// \returns A pair of two vectors. The first is a concatenation of all input strings
/// (without any separating null characters). The second is the list of lengths of these strings.
std::pair<std::vector<char>, std::vector<size_t>> encodeStrings(
        const std::vector<std::string> &strings) {
    std::pair<std::vector<char>, std::vector<size_t>> result;
    std::vector<char> &charArray = result.first;
    std::vector<size_t> &lengths = result.second;

    size_t totalLength = 0;
    lengths.reserve(strings.size());
    for (const std::string &s : strings) {
        lengths.push_back(s.size());
        totalLength += s.size();
    }

    charArray.reserve(totalLength);
    for (const std::string &s : strings) {
        charArray.insert(charArray.end(), s.begin(), s.end());
    }

    return result;
}

/// \brief Split a character array into multiple strings.
///
/// \param charArray
///   A character array storing a number of concatenated strings (without separating null
///   characters).
///
/// \param lengths
///  The list of lengths of the strings stored in \p charArray.
///
/// \returns A vector of strings extracted from \p charArray.
std::vector<std::string> decodeStrings(const std::vector<char> &charArray,
                                       const std::vector<size_t> &lengths) {
    std::vector<std::string> strings;
    strings.reserve(lengths.size());

    std::vector<char>::const_iterator nextStringBegin = charArray.begin();
    for (size_t length : lengths) {
        strings.emplace_back(nextStringBegin, nextStringBegin + length);
        nextStringBegin += length;
    }

    return strings;
}

}  // namespace

namespace oops {
namespace mpi {

// ------------------------------------------------------------------------------------------------

const eckit::mpi::Comm & clone(const eckit::mpi::Comm & comm) {
  return eckit::mpi::comm(comm.name().c_str());
}

// ------------------------------------------------------------------------------------------------

const eckit::mpi::Comm & world() {
  return eckit::mpi::comm("world");
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

void allGatherv(const eckit::mpi::Comm & comm, std::vector<util::DateTime> &x) {
    size_t globalSize = x.size();
    comm.allReduceInPlace(globalSize, eckit::mpi::sum());
    std::vector<util::DateTime> globalX(globalSize);
    oops::mpi::allGathervUsingSerialize(comm, x.begin(), x.end(), globalX.begin());
    x = std::move(globalX);
}

// ------------------------------------------------------------------------------------------------

void allGatherv(const eckit::mpi::Comm & comm, std::vector<std::string> &x) {
    std::pair<std::vector<char>, std::vector<size_t>> encodedX = encodeStrings(x);

    // Gather all character arrays
    eckit::mpi::Buffer<char> charBuffer(comm.size());
    comm.allGatherv(encodedX.first.begin(), encodedX.first.end(), charBuffer);

    // Gather all string lengths
    eckit::mpi::Buffer<size_t> lengthBuffer(comm.size());
    comm.allGatherv(encodedX.second.begin(), encodedX.second.end(), lengthBuffer);

    // Free memory
    encodedX = {};

    x = decodeStrings(charBuffer.buffer, lengthBuffer.buffer);
}

// ------------------------------------------------------------------------------------------------

void exclusiveScan(const eckit::mpi::Comm &comm, size_t &x) {
  // Could be done with MPI_Exscan, but there's no wrapper for it in eckit::mpi.

  std::vector<size_t> xs(comm.size());
  comm.allGather(x, xs.begin(), xs.end());
  x = std::accumulate(xs.begin(), xs.begin() + comm.rank(), 0);
}

// ------------------------------------------------------------------------------------------------
void broadcastBool(const eckit::mpi::Comm & comm, bool & boolVar, const int root) {
    // Send bool as int since eckit MPI broadcast doesn't accept bool
    int tempInt;
    if (comm.rank() == root) {
        tempInt = static_cast<int>(boolVar);
        comm.broadcast(tempInt, root);
    } else {
        comm.broadcast(tempInt, root);
        boolVar = static_cast<bool>(tempInt);
    }
}

void broadcastString(const eckit::mpi::Comm & comm, std::string & stringVar, const int root) {
    std::vector<char> buffer;
    if (comm.rank() == root) {
        // Send string as vector of char since eckit MPI broadcast doesn't accept string
        buffer.resize(stringVar.size() + 1);  // allow for trailing NULL
        strncpy(buffer.data(), stringVar.data(), stringVar.size());
        buffer[stringVar.size()] = '\0';      // make sure there is a trailing NULL
        broadcastVector<char>(comm, buffer, root);
    } else {
        broadcastVector<char>(comm, buffer, root);
        stringVar = buffer.data();
    }
}

// ------------------------------------------------------------------------------------------------

// MPI tag values for send/receive utilities
static const int msgIsSize = 1;
static const int msgIsData = 2;

void sendString(const eckit::mpi::Comm & comm, const std::string & stringVar,
                                               const int toRank)  {
    // First send the string length, then send the string
    int stringSize = stringVar.size();
    std::vector<char> buffer(stringSize + 1);    // Allow for trailing NULL
    strncpy(buffer.data(), stringVar.data(), stringSize);
    buffer[stringSize] = '\0';
    comm.send<int>(&stringSize, 1, toRank, msgIsSize);
    comm.send<char>(buffer.data(), buffer.size(), toRank, msgIsData);
}

void receiveString(const eckit::mpi::Comm & comm, std::string & stringVar,
                                                  const int fromRank) {
    // First receive the string length, then receive the string
    // The string will be coming in a vector of char
    int stringSize;
    comm.receive<int>(&stringSize, 1, fromRank, msgIsSize);
    // Can't assign directly to a string.data() pointer (which is const);
    std::vector<char> buffer(stringSize + 1);
    comm.receive<char>(buffer.data(), buffer.size(), fromRank, msgIsData);
    stringVar = buffer.data();
}

}  // namespace mpi
}  // namespace oops
