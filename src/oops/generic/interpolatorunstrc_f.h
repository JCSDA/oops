/*
 * (C) Copyright 2020- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include "atlas/field.h"

// Forward declarations
namespace eckit {
  class Configuration;
  namespace mpi {
    class Comm;
  }
}

namespace oops {
extern "C" {
  void unstrc_create_f90(int &, const eckit::mpi::Comm *,
                       const atlas::field::FieldImpl *,
                       const atlas::field::FieldImpl *,
                       const eckit::Configuration &);
  void unstrc_write_f90(const int &, const eckit::Configuration &);
  void unstrc_delete_f90(const int &);
  void unstrc_apply_f90(const int &, const atlas::field::FieldImpl *,
                       atlas::field::FieldImpl *);
  void unstrc_apply_ad_f90(const int &, const atlas::field::FieldImpl *,
                       atlas::field::FieldImpl *);
}
}  // namespace oops
