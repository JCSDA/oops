/*
 * (C) Copyright 2020- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "eckit/config/Configuration.h"

#include "oops/base/InterpolatorBase.h"
#include "oops/mpi/mpi.h"
#include "oops/util/ObjectCounter.h"

namespace oops {

// -----------------------------------------------------------------------------

/*! \brief Interface for Unstructured interpolation
 *
 */

class InterpolatorUnstructured : public InterpolatorBase,
                                 private util::ObjectCounter<InterpolatorUnstructured> {
 public:
  static const std::string classname() {return "oops::InterpolatorUnstructured";}

  InterpolatorUnstructured(const eckit::Configuration &, const atlas::FunctionSpace &,
                           const atlas::FunctionSpace &,
                           const atlas::field::FieldSetImpl * = nullptr,
                           const eckit::mpi::Comm & = oops::mpi::world());
  ~InterpolatorUnstructured();

  void apply(const atlas::Field &, atlas::Field &) override;
  void apply(const atlas::FieldSet &, atlas::FieldSet &) override;

  void apply_ad(const atlas::Field &, atlas::Field &);
  void apply_ad(const atlas::FieldSet &, atlas::FieldSet &) override;

  int write(const eckit::Configuration &) override;

 private:
  int keyUnstructuredInterpolator_;
  const atlas::FunctionSpace *in_fspace_;
  const atlas::FunctionSpace *out_fspace_;
  void print(std::ostream &) const override;
};

// -----------------------------------------------------------------------------

}  // namespace oops
