/*
 * (C) Crown Copyright 2023 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <netcdf.h>
#include <string>
#include <type_traits>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "atlas/array.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/missingValues.h"

namespace eckit {
namespace mpi {
class Comm;
}  // namespace mpi
}  // namespace eckit

// Assumes that ArrayView shape and value is the same on all PEs
//
#define ERR1(e) {::util::abor1_cpp(nc_strerror(e));}

namespace util {

/// \brief atlas array gather for 1 dimensional ArrayView
///        and sum the contributions from all the pes on the root PE
template<typename Value>
void gatherSum(const eckit::mpi::Comm & comm,
               const std::size_t & root,
               atlas::array::ArrayView<Value, 1> & arrayInOut);

/// \brief atlas array gather for 2 dimensional ArrayView
///        and sum the contributions from all the pes on the root PE
template<typename Value>
void gatherSum(const eckit::mpi::Comm & comm,
               const std::size_t & root,
               atlas::array::ArrayView<Value, 2> & arrayInOut);

/// \brief atlas array gather for 3 dimensional ArrayView
///        and sum the contributions from all the pes on the root PE
template<typename Value>
void gatherSum(const eckit::mpi::Comm & comm,
               const std::size_t & root,
               atlas::array::ArrayView<Value, 3> & arrayInOut);

/// \brief atlas array scatter for 1 dimensional ArrayView
///        from root PE to all pes.
template<typename Value>
void scatter(const eckit::mpi::Comm & comm,
             const std::size_t & root,
             atlas::array::ArrayView<Value, 1> & arrayInOut);

/// \brief atlas array scatter for 2 dimensional ArrayView
///        from root PE to all pes.
template<typename Value>
void scatter(const eckit::mpi::Comm & comm,
             const std::size_t & root,
             atlas::array::ArrayView<Value, 2> & arrayInOut);

/// \brief atlas array scatter for 3 dimensional ArrayView
///        from root PE to all pes.
template<typename Value>
void scatter(const eckit::mpi::Comm & comm,
             const std::size_t & root,
             atlas::array::ArrayView<Value, 3> & arrayInOut);

/// \brief - creates the header in the written NetCDF file
///        - assumes that we have the:
///        - netcdf file path (ncfilepath)
///        - dimension names (dimNames) for all dimensions in the file to be written
///        - dimension sizes for each of the dimension names. (Assumed to be in the same
///          order as dimNames
///        - the variableNames of the fields in the file
///        - the dimension names for each variable in C++ ordering
// assuming the netcdf_ids on input has zero size.
void atlasArrayWriteHeader(
    const std::string & ncfilepath,
    const std::vector<std::string> & dimNames,
    const std::vector<atlas::idx_t> & dimSizes,
    const std::vector<std::string> & variableNames,
    const std::vector<std::vector<std::string>> & dimNamesForEveryVar,
    std::vector<int> & netcdfGeneralIDs,
    std::vector<int> & netcdfDimIDs,
    std::vector<int> & netcdfVarIDs,
    std::vector<std::vector<int>> & netcdf_dim_varIDs);

/// \brief - interrogates the header of a written NetCDF file
///        - assumes that we have the:
///        - netcdf file path (ncfilepath)
void atlasArrayInquire(
    const std::string & ncfilepath,
    std::vector<std::string> & dimNames,
    std::vector<atlas::idx_t> & dimSizes,
    std::vector<std::string> & variableNames,
    std::vector<std::vector<std::string>> & dimNamesForEveryVar,
    std::vector<int> & netcdfGeneralIDs,
    std::vector<int> & netcdfDimIDs,
    std::vector<int> & netcdfVarIDs,
    std::vector<std::vector<int>> & netcdf_dim_varIDs);

/// \brief - reads variable data that is 1 dimensional and of double type
void atlasArrayReadData(
    const std::vector<int> & netcdfGeneralIDs,
    const std::vector<atlas::idx_t> & dimSizes,
    const int & varID,
    atlas::array::ArrayView<double, 1> & arrayInOut);

/// \brief - reads variable data that is 2 dimensional and of double type
void atlasArrayReadData(
    const std::vector<int> & netcdfGeneralIDs,
    const std::vector<atlas::idx_t> & dimSizes,
    const int & varID,
    atlas::array::ArrayView<double, 2> & arrayInOut);

/// \brief - reads variable data that is 3 dimensional and of double type
void atlasArrayReadData(
    const std::vector<int> & netcdfGeneralIDs,
    const std::vector<atlas::idx_t> & dimSizes,
    const int & varID,
    atlas::array::ArrayView<double, 3> & arrayInOut);

/// \brief - reads variable data that is 1 dimensional and of integer type
void atlasArrayReadData(
    const std::vector<int> & netcdfGeneralIDs,
    const std::vector<atlas::idx_t> & dimSizes,
    const int & varID,
    atlas::array::ArrayView<int, 1> & arrayInOut);

/// \brief - writes variable data that is 1 dimensional and of double type
void atlasArrayWriteData(
    const std::vector<int> & netcdfGeneralIDs,
    const int & varID,
    atlas::array::ArrayView<const double, 1> & arrayIn);

/// \brief - writes variable data that is 2 dimensional and of double type
void atlasArrayWriteData(
    const std::vector<int> & netcdfGeneralIDs,
    const int & varID,
    atlas::array::ArrayView<const double, 2> & arrayIn);

/// \brief - writes variable data that is 3 dimensional and of double type
void atlasArrayWriteData(
    const std::vector<int> & netcdfGeneralIDs,
    const int & varID,
    atlas::array::ArrayView<const double, 3> & arrayIn);

/// \brief - writes variable data that is 3 dimensional and of integer type
void atlasArrayWriteData(
    const std::vector<int> & netcdfGeneralIDs,
    const int & varID,
    atlas::array::ArrayView<const int, 1> & arrayIn);

//-----------------------------------------------------------------------------
// Implementation
//-----------------------------------------------------------------------------
using atlas::idx_t;


template<typename Value>
void gatherSum(const eckit::mpi::Comm & comm,
               const std::size_t & root,
               atlas::array::ArrayView<Value, 1> & arrayInOut) {
  std::vector<Value> recv(arrayInOut.size() * comm.size(),
                          static_cast<Value>(0));
  std::vector<Value> send(arrayInOut.size(), static_cast<Value>(0));

  size_t n = 0;
  for (idx_t t = 0; t < arrayInOut.shape()[0]; ++t) {
    send[n] = arrayInOut(t);
    ++n;
  }

  comm.gather(send, recv, root);

  if (comm.rank() == root) {
    n = 0;
    size_t stride = arrayInOut.size();
    for (idx_t t = 0; t < arrayInOut.shape()[0]; ++t) {
      arrayInOut(t) = static_cast<Value>(0);
      for (size_t pe = 0; pe < comm.size(); ++pe) {
          arrayInOut(t) += recv[n + pe * stride];
      }
      ++n;
    }
  }
}


template<typename Value>
void gatherSum(const eckit::mpi::Comm & comm,
               const std::size_t & root,
               atlas::array::ArrayView<Value, 2> & arrayInOut) {
  std::vector<Value> recv(arrayInOut.size() * comm.size(),
                          static_cast<Value>(0));
  std::vector<Value> send(arrayInOut.size(), static_cast<Value>(0));

  size_t n = 0;
  for (idx_t t = 0; t < arrayInOut.shape()[0]; ++t) {
    for (idx_t t1 = 0; t1 < arrayInOut.shape()[1]; ++t1) {
      send[n] = arrayInOut(t, t1);
      ++n;
    }
  }

  comm.gather(send, recv, root);

  if (comm.rank() == root) {
    n = 0;
    size_t stride = arrayInOut.size();
    for (idx_t t = 0; t < arrayInOut.shape()[0]; ++t) {
      for (idx_t t1 = 0; t1 < arrayInOut.shape()[1]; ++t1) {
        arrayInOut(t, t1) = static_cast<Value>(0);
        for (size_t pe = 0; pe < comm.size(); ++pe) {
          arrayInOut(t, t1) += recv[n + pe * stride];
        }
        ++n;
      }
    }
  }
}


template<typename Value>
void gatherSum(const eckit::mpi::Comm & comm,
               const std::size_t & root,
               atlas::array::ArrayView<Value, 3> & arrayInOut) {
  std::vector<Value> recv(arrayInOut.size() * comm.size(),
                          static_cast<Value>(0));
  std::vector<Value> send(arrayInOut.size(), static_cast<Value>(0));

  size_t n = 0;
  for (idx_t t = 0; t < arrayInOut.shape()[0]; ++t) {
    for (idx_t t1 = 0; t1 < arrayInOut.shape()[1]; ++t1) {
      for (idx_t t2 = 0; t2 < arrayInOut.shape()[2]; ++t2) {
        send[n] = arrayInOut(t, t1, t2);
        ++n;
      }
    }
  }

  comm.gather(send, recv, root);

  if (comm.rank() == root) {
    n = 0;
    size_t stride = arrayInOut.size();
    for (idx_t t = 0; t < arrayInOut.shape()[0]; ++t) {
      for (idx_t t1 = 0; t1 < arrayInOut.shape()[1]; ++t1) {
        for (idx_t t2 = 0; t2 < arrayInOut.shape()[2]; ++t2) {
          arrayInOut(t, t1, t2) = static_cast<Value>(0);
          for (size_t pe = 0; pe < comm.size(); ++pe) {
            arrayInOut(t, t1, t2) += recv[n + pe * stride];
          }
          ++n;
        }
      }
    }
  }
}


template<typename Value>
void scatter(
    const eckit::mpi::Comm & comm,
    const std::size_t & root,
    atlas::array::ArrayView<Value, 1> & arrayInOut) {
  std::vector<Value> recv(arrayInOut.size(), static_cast<Value>(0));
  std::vector<Value> send(arrayInOut.size() * comm.size(), static_cast<Value>(0));

  if (comm.rank() == root) {
    for (idx_t pe = 0; pe < static_cast<idx_t>(comm.size()); ++pe) {
      size_t n = 0;
      for (idx_t t = 0; t < arrayInOut.shape()[0]; ++t) {
        send[arrayInOut.size() * pe + n] = arrayInOut(t);
        ++n;
      }
    }
  }

  comm.scatter(send, recv, root);

  size_t n = 0;
  for (idx_t t = 0; t < arrayInOut.shape()[0]; ++t) {
    arrayInOut(t) = recv[n];
    ++n;
  }
}


template<typename Value>
void scatter(
    const eckit::mpi::Comm & comm,
    const std::size_t & root,
    atlas::array::ArrayView<Value, 2> & arrayInOut) {
  std::vector<Value> recv(arrayInOut.size(), static_cast<Value>(0));
  std::vector<Value> send(arrayInOut.size() * comm.size(), static_cast<Value>(0));

  if (comm.rank() == root) {
    for (idx_t pe = 0; pe < static_cast<idx_t>(comm.size()); ++pe) {
      size_t n = 0;
      for (idx_t t = 0; t < arrayInOut.shape()[0]; ++t) {
        for (idx_t t1 = 0; t1 < arrayInOut.shape()[1]; ++t1) {
          send[arrayInOut.size() * pe + n] = arrayInOut(t, t1);
          ++n;
        }
      }
    }
  }

  comm.scatter(send, recv, root);

  size_t n = 0;
  for (idx_t t = 0; t < arrayInOut.shape()[0]; ++t) {
    for (idx_t t1 = 0; t1 < arrayInOut.shape()[1]; ++t1) {
      arrayInOut(t, t1) = recv[n];
      ++n;
    }
  }
}


template<typename Value>
void scatter(
    const eckit::mpi::Comm & comm,
    const std::size_t & root,
    atlas::array::ArrayView<Value, 3> & arrayInOut) {
  std::vector<Value> recv(arrayInOut.size(), static_cast<Value>(0));
  std::vector<Value> send(arrayInOut.size() * comm.size(), static_cast<Value>(0));

  if (comm.rank() == root) {
    for (idx_t pe = 0; pe < static_cast<idx_t>(comm.size()); ++pe) {
      size_t n = 0;
      for (idx_t t = 0; t < arrayInOut.shape()[0]; ++t) {
        for (idx_t t1 = 0; t1 < arrayInOut.shape()[1]; ++t1) {
          for (idx_t t2 = 0; t2 < arrayInOut.shape()[2]; ++t2) {
            send[arrayInOut.size() * pe + n] = arrayInOut(t, t1, t2);
            ++n;
          }
        }
      }
    }
  }

  comm.scatter(send, recv, root);

  size_t n = 0;
  for (idx_t t = 0; t < arrayInOut.shape()[0]; ++t) {
    for (idx_t t1 = 0; t1 < arrayInOut.shape()[1]; ++t1) {
      for (idx_t t2 = 0; t2 < arrayInOut.shape()[2]; ++t2) {
        arrayInOut(t, t1, t2) = recv[n];
        ++n;
      }
    }
  }
}

}  // namespace util
