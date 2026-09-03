/*
 * (C) Copyright 2022- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "oops/base/GeometryData.h"

#include "eckit/mpi/Comm.h"

namespace oops {

// -----------------------------------------------------------------------------

GeometryData::GeometryData(const atlas::FunctionSpace & fspace, const atlas::FieldSet & fset,
                           const bool topdown, const eckit::mpi::Comm & comm):
  fspace_(fspace), fset_(fset), comm_(comm), topdown_(topdown)
{
  // Set default communicator name
  eckit::mpi::setCommDefault(comm_.name().c_str());
}

// -----------------------------------------------------------------------------

}  // namespace oops
