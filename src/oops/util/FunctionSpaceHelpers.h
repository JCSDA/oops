/*
 * (C) Copyright 2024 UCAR
 * (C) Crown Copyright 2024 Met Office
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <utility>
#include <vector>

#include "atlas/functionspace.h"

#include "oops/util/Logger.h"

namespace atlas {
namespace grid { class Distribution; }
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

/// \brief procedure to call a functor for a given concrete implementation
///        of a function space type
template<typename Functor>
void executeFunc(const atlas::FunctionSpace & fspace, const Functor & functor) {
  if (atlas::functionspace::NodeColumns(fspace)) {
    functor(atlas::functionspace::CubedSphereNodeColumns(fspace));
  } else if (atlas::functionspace::StructuredColumns(fspace)) {
    functor(atlas::functionspace::StructuredColumns(fspace));
  } else {
    oops::Log::error() << "ERROR - a functor call failed "
                          "(function space type not allowed)" << std::endl;
    throw std::runtime_error("a functor call failed");
  }
}

atlas::idx_t getSizeOwned(const atlas::FunctionSpace & fspace);

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

// Define a distribution and a mesh from a custom partition
void setupStructuredMeshWithCustomPartition(const eckit::mpi::Comm &,
    const atlas::Grid &,
    const std::vector<int> &,
    atlas::grid::Distribution &,
    atlas::Mesh &);

// -----------------------------------------------------------------------------

}  // namespace util
