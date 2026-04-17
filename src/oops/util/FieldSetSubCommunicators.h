/*
 * (C) Crown Copyright 2024, Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field.h"
#include "atlas/functionspace.h"

#include "eckit/mpi/Comm.h"

#include "oops/util/redistribution/CommRedistribution.h"

namespace util {

void redistributeToSubcommunicator(const CommRedistribution & redist,
                                   const atlas::FieldSet & fsetIn,
                                   atlas::FieldSet & fsetOut,
                                   const atlas::FunctionSpace & fspaceOut);

/// Copy fieldSet from communicator to sub-communicators
void redistributeToSubcommunicator(const atlas::FieldSet &,
                                   atlas::FieldSet &,
                                   const eckit::mpi::Comm &,
                                   const eckit::mpi::Comm &,
                                   const atlas::FunctionSpace &,
                                   const atlas::FunctionSpace &);


void gatherAndSumFromSubcommunicator(const CommRedistribution & redist,
                                     const atlas::FieldSet & fsetIn,
                                     atlas::FieldSet & fsetOut,
                                     const eckit::mpi::Comm & comm,
                                     const atlas::FunctionSpace & fspaceIn,
                                     const atlas::FunctionSpace & fspaceOut);

/// Gather and sum fieldSets from sub-communicators to larger communicator
void gatherAndSumFromSubcommunicator(const atlas::FieldSet &,
                                     atlas::FieldSet &,
                                     const eckit::mpi::Comm &,
                                     const eckit::mpi::Comm &,
                                     const atlas::FunctionSpace &,
                                     const atlas::FunctionSpace &);
}  // namespace util
