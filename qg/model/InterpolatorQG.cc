/*
 * (C) Copyright 2022-2022 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "model/InterpolatorQG.h"

#include <ostream>
#include <vector>

#include "model/GeometryQG.h"
#include "model/IncrementQG.h"
#include "model/QgFortran.h"
#include "model/StateQG.h"

namespace qg {

// -----------------------------------------------------------------------------

InterpolatorQG::InterpolatorQG(const eckit::Configuration &,
                               const GeometryQG & grid, const std::vector<double> & locs)
  : nlevs_(grid.levels()), nlocs_(locs.size() / 2), locs_(locs)
{
  ASSERT(locs.size() % 2 == 0);
}

// -----------------------------------------------------------------------------

InterpolatorQG::~InterpolatorQG() {}

// -----------------------------------------------------------------------------

void InterpolatorQG::apply(const oops::Variables & vars, const StateQG & xx,
                           std::vector<double> & values) const {
  const size_t nvals = vars.size() * nlevs_ * nlocs_;
  values.resize(nvals);
  qg_fields_getvals_f90(xx.fields().toFortran(), vars, nlocs_, locs_[0], nvals, values[0]);
}

// -----------------------------------------------------------------------------

void InterpolatorQG::apply(const oops::Variables & vars, const IncrementQG & dx,
                           std::vector<double> & values) const {
  const size_t nvals = vars.size() * nlevs_ * nlocs_;
  values.resize(nvals);
  qg_fields_getvals_f90(dx.fields().toFortran(), vars, nlocs_, locs_[0], nvals, values[0]);
}

// -----------------------------------------------------------------------------

void InterpolatorQG::applyAD(const oops::Variables & vars, IncrementQG & dx,
                             const std::vector<double> & values) const {
  const size_t nvals = vars.size() * nlevs_ * nlocs_;
  ASSERT(values.size() == nvals);
  qg_fields_getvalsad_f90(dx.fields().toFortran(), vars, nlocs_, locs_[0], nvals, values[0]);
}

// -----------------------------------------------------------------------------

void InterpolatorQG::print(std::ostream & os) const {
  os << "InterpolatorQG";
}

// -----------------------------------------------------------------------------

}  // namespace qg

