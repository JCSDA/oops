/*
 * (C) Copyright 2022- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"

namespace eckit {
namespace mpi {
class Comm;
}  // namespace mpi
}  // namespace eckit

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Lightweight wrapper over a model geometry as represented via Atlas.
class GeometryData {
 public:
  GeometryData(const atlas::FunctionSpace &, const atlas::FieldSet &,
               const bool, const eckit::mpi::Comm &);

  ~GeometryData() = default;

  GeometryData(const GeometryData &) = delete;
  GeometryData & operator=(const GeometryData &) = delete;

  const atlas::FunctionSpace & functionSpace() const {return fspace_;}
  const atlas::FieldSet & fieldSet() const {return fset_;}
  const eckit::mpi::Comm & comm() const {return comm_;}

  bool has(const std::string & name) const {return fset_.has(name);}
  const atlas::Field & getField(const std::string & name) const {return fset_.field(name);}
  bool levelsAreTopDown() const {return topdown_;}

 private:
  atlas::FunctionSpace fspace_;
  atlas::FieldSet fset_;
  const eckit::mpi::Comm & comm_;
  bool topdown_;
};

// -----------------------------------------------------------------------------

}  // namespace oops
