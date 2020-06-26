/*
 * (C) Copyright 2020- UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <ostream>
#include <string>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace.h"
#include "atlas/interpolation.h"
#include "atlas/option.h"
#include "eckit/config/Configuration.h"

#include "oops/generic/InterpolatorAtlas.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"

using atlas::array::make_datatype;
using atlas::array::make_shape;
using atlas::option::name;
using atlas::option::levels;

namespace oops {

static InterpolatorMaker<InterpolatorAtlas> makerId_("atlas");

// -----------------------------------------------------------------------------
InterpolatorAtlas::InterpolatorAtlas(const eckit::Configuration & config,
                   const atlas::FunctionSpace & grid1,
                   const atlas::FunctionSpace & grid2,
                   const atlas::field::FieldSetImpl * masks) {
    // for now, let's use the finite element interpolation, which is available
    // for unstructured grids.  We may wish to implement more options in the
    // future.  Also, atlas has not yet implemented masks so ignore that
    // argument.
      interpolator_.reset(new atlas::Interpolation(atlas::option::type( "finite-element" ),
                                                   grid1, grid2));
}

// -----------------------------------------------------------------------------
void InterpolatorAtlas::apply(atlas::Field const & infield,
                              atlas::Field & outfield) {
      interpolator_->execute(infield, outfield);
}

// -----------------------------------------------------------------------------
void InterpolatorAtlas::apply(atlas::FieldSet const & infields,
                              atlas::FieldSet & outfields) {
  // Allocate space for the output fields if the caller has not already done so
  for (size_t ifield = 0; ifield < static_cast<size_t>(infields.size()); ++ifield) {
    std::string fname = infields.field(ifield).name();
    if (!outfields.has_field(fname)) {
      oops::Log::info() << "Allocating output fields for Atlas Interpolation" << std::endl;

      atlas::Field outfield =
        interpolator_->target().createField<double>(name(fname) |
                                levels(infields.field(ifield).levels()));
      outfields.add(outfield);
    }
  }

  // apply interpolation
  interpolator_->execute(infields, outfields);
}

// -----------------------------------------------------------------------------
void InterpolatorAtlas::print(std::ostream & os) const {
  os << " InterpolatorAtlas: print not implemented yet.";
}
// -----------------------------------------------------------------------------
}  // namespace oops
