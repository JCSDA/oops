/*
 * (C) Copyright 2024 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

namespace atlas {
namespace grid { class Partitioner; }
class Grid;
class Mesh;
class FunctionSpace;
class FieldSet;
}  // namespace atlas

namespace eckit {
namespace mpi { class Comm; }
class Configuration;
}  // namespace eckit

namespace util {

// -----------------------------------------------------------------------------

// Parses config and sets up an atlas::FunctionSpace and associated geometric data
void setupFunctionSpace(const eckit::mpi::Comm & comm,
    const eckit::Configuration & config,
    atlas::Grid & grid,
    atlas::grid::Partitioner & partitioner,
    atlas::Mesh & mesh,
    atlas::FunctionSpace & functionspace,
    atlas::FieldSet & fieldset);

// -----------------------------------------------------------------------------

}  // namespace util
