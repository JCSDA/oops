/*
 * (C) Copyright 2020- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/generic/interpolatorunstrc_f.h"
#include "oops/generic/InterpolatorUnstructured.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

namespace oops {

static InterpolatorMaker<InterpolatorUnstructured> makerUNSTR_("unstructured");

// -----------------------------------------------------------------------------
InterpolatorUnstructured::InterpolatorUnstructured(const eckit::Configuration & config,
                                                   const atlas::FunctionSpace & fspace1,
                                                   const atlas::FunctionSpace & fspace2,
                                                   const atlas::field::FieldSetImpl * masks,
                                                   const eckit::mpi::Comm & comm)
  : in_fspace_(&fspace1), out_fspace_(&fspace2)
{
  // mask have not yet been implemented for unstructured interpolation
  unstrc_create_f90(keyUnstructuredInterpolator_, &comm,
                    fspace1.lonlat().get(), fspace2.lonlat().get(), config);
}

// -----------------------------------------------------------------------------
int InterpolatorUnstructured::write(const eckit::Configuration & config) {
  unstrc_write_f90(keyUnstructuredInterpolator_, config);
  return 0;
}

// -----------------------------------------------------------------------------
InterpolatorUnstructured::~InterpolatorUnstructured() {
  unstrc_delete_f90(keyUnstructuredInterpolator_);
}
// -----------------------------------------------------------------------------
void InterpolatorUnstructured::apply(const atlas::Field & infield, atlas::Field & outfield) {
  unstrc_apply_f90(keyUnstructuredInterpolator_, infield.get(), outfield.get());
}

// -----------------------------------------------------------------------------
void InterpolatorUnstructured::apply_ad(const atlas::Field & field_grid2,
                                        atlas::Field & field_grid1) {
  unstrc_apply_ad_f90(keyUnstructuredInterpolator_, field_grid2.get(), field_grid1.get());
}

// -----------------------------------------------------------------------------
void InterpolatorUnstructured::apply(const atlas::FieldSet & infields,
                                     atlas::FieldSet & outfields) {
  // Allocate space for the output fields if the caller has not already done so
  for (int ifield = 0; ifield < infields.size(); ++ifield) {
    std::string fname = infields.field(ifield).name();
    if (!outfields.has_field(fname)) {
      atlas::Field outfield = out_fspace_->createField<double>(atlas::option::name(fname) |
                              atlas::option::levels(infields.field(ifield).levels()));
      outfields.add(outfield);
    }
    this->apply(infields.field(fname), outfields.field(fname));
  }
}

// -----------------------------------------------------------------------------
void InterpolatorUnstructured::apply_ad(const atlas::FieldSet & fields_grid2,
                                        atlas::FieldSet & fields_grid1) {
  // Allocate space for the output fields if the caller has not already done so
  for (int ifield = 0; ifield < fields_grid2.size(); ++ifield) {
    std::string fname = fields_grid2.field(ifield).name();
    if (!fields_grid1.has_field(fname)) {
      atlas::Field field1 = in_fspace_->createField<double>(atlas::option::name(fname) |
                            atlas::option::levels(fields_grid2.field(ifield).levels()));
      fields_grid1.add(field1);
    }
    this->apply_ad(fields_grid2.field(fname), fields_grid1.field(fname));
  }
}

// -----------------------------------------------------------------------------
void InterpolatorUnstructured::print(std::ostream & os) const {
  os << " InterpolatorUnstructured: print not implemented yet.";
}
// -----------------------------------------------------------------------------
}  // namespace oops
