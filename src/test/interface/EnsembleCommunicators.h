/*
 * (C) Copyright 2026 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef TEST_INTERFACE_ENSEMBLECOMMUNICATORS_H_
#define TEST_INTERFACE_ENSEMBLECOMMUNICATORS_H_

#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"
#include "oops/mpi/mpi.h"

namespace test {

// -----------------------------------------------------------------------------
/// Communicators of a test's tasks reinterpreted as an ensemble
struct EnsembleCommunicators {
  const eckit::mpi::Comm & member;  ///< tasks holding my ensemble member
  const eckit::mpi::Comm & ens;     ///< tasks holding my patch of every member
};

// -----------------------------------------------------------------------------
/*! \brief Reinterprets the test's T tasks as N ensemble members of T/N tasks each
 *
 * \details N is set by the option `allreduce members` in \p testConf (default T, i.e. each
 * member on 1 task), and must divide T. The communicators mirror LocalEnsembleDA: members are
 * contiguous blocks of tasks, each with its own communicator (<prefix>comm_member_<i>), and the
 * tasks holding the same patch of every member share an ensemble communicator
 * (<prefix>patch_member_<j>). Communicators that already exist (e.g. created by an earlier call
 * with the same prefix) are reused.
 */
inline EnsembleCommunicators splitEnsembleCommunicators(const eckit::Configuration & testConf,
                                                        const std::string & prefix) {
  const eckit::mpi::Comm & world = oops::mpi::world();
  const int ntasks = static_cast<int>(world.size());
  const int nmembers = testConf.getInt("allreduce members", ntasks);
  if (nmembers <= 0 || ntasks % nmembers != 0) {
    throw eckit::BadValue("allreduce members (" + std::to_string(nmembers)
                          + ") must divide the number of tasks ("
                          + std::to_string(ntasks) + ")", Here());
  }
  const int tasksPerMember = ntasks / nmembers;
  const int myMember = static_cast<int>(world.rank()) / tasksPerMember;
  const int myPatch = static_cast<int>(world.rank()) % tasksPerMember;
  const std::string memberName = prefix + "comm_member_" + std::to_string(myMember);
  const std::string ensName = prefix + "patch_member_" + std::to_string(myPatch);
  const eckit::mpi::Comm & commMember = eckit::mpi::hasComm(memberName)
                                        ? eckit::mpi::comm(memberName)
                                        : world.split(myMember, memberName);
  const eckit::mpi::Comm & commEns = eckit::mpi::hasComm(ensName)
                                     ? eckit::mpi::comm(ensName)
                                     : world.split(myPatch, ensName);
  ASSERT(static_cast<int>(commMember.size()) == tasksPerMember);
  ASSERT(static_cast<int>(commEns.size()) == nmembers);
  return EnsembleCommunicators{commMember, commEns};
}

// -----------------------------------------------------------------------------

}  // namespace test

#endif  // TEST_INTERFACE_ENSEMBLECOMMUNICATORS_H_
