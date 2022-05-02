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

InterpolatorQG::InterpolatorQG(const eckit::Configuration &, const GeometryQG & grid,
                               const std::vector<double> & lats, const std::vector<double> & lons)
  : nlevs_(grid.levels()), nlocs_(lats.size()), locs_(2 * nlocs_)
{
  ASSERT(lats.size() == lons.size());
  for (size_t jj = 0; jj < nlocs_; ++jj) {
    locs_[2 * jj] = lats[jj];
    locs_[2 * jj + 1] = lons[jj];
  }
}

// -----------------------------------------------------------------------------

InterpolatorQG::~InterpolatorQG() {}

// -----------------------------------------------------------------------------

void InterpolatorQG::apply(const oops::Variables & vars, const StateQG & xx,
                           const std::vector<bool> & mask,
                           std::vector<double> & values) const {
  this->apply(vars, xx.fields(), mask, values);
}

// -----------------------------------------------------------------------------

void InterpolatorQG::apply(const oops::Variables & vars, const IncrementQG & dx,
                           const std::vector<bool> & mask,
                           std::vector<double> & values) const {
  this->apply(vars, dx.fields(), mask, values);
}

// -----------------------------------------------------------------------------

void InterpolatorQG::apply(const oops::Variables & vars, const FieldsQG & flds,
                           const std::vector<bool> & mask,
                           std::vector<double> & values) const {
  ASSERT(mask.size() == nlocs_);
  ASSERT(values.size() == vars.size() * nlevs_ * nlocs_);

  size_t nout = 0;
  std::vector<double> locs;
  for (size_t jj = 0; jj < nlocs_; ++jj) {
    if (mask[jj]) {
      locs.push_back(locs_[2*jj]);
      locs.push_back(locs_[2*jj+1]);
      ++nout;
    }
  }
  const size_t nvals = vars.size() * nlevs_ * nout;
  std::vector<double> tmp(nvals);

  qg_fields_getvals_f90(flds.toFortran(), vars, nout, locs[0], nvals, tmp[0]);

  size_t itmp = 0;
  size_t ival = 0;
  for (size_t jv = 0; jv < vars.size(); ++jv) {
    for (size_t jj = 0; jj < nlocs_; ++jj) {
      if (mask[jj]) {
        for (size_t jl = 0; jl < nlevs_; ++jl) {
          values[ival + jl] = tmp[itmp + jl];
        }
        itmp += nlevs_;
      }
      ival += nlevs_;
    }
  }
}

// -----------------------------------------------------------------------------

void InterpolatorQG::applyAD(const oops::Variables & vars, IncrementQG & dx,
                             const std::vector<bool> & mask,
                             const std::vector<double> & values) const {
  ASSERT(mask.size() == nlocs_);
  ASSERT(values.size() == vars.size() * nlevs_ * nlocs_);

  size_t nout = 0;
  std::vector<double> locs;
  for (size_t jj = 0; jj < nlocs_; ++jj) {
    if (mask[jj]) {
      locs.push_back(locs_[2*jj]);
      locs.push_back(locs_[2*jj+1]);
      ++nout;
    }
  }
  const size_t nvals = vars.size() * nlevs_ * nout;
  std::vector<double> tmp(nvals);

  size_t itmp = 0;
  size_t ival = 0;
  for (size_t jv = 0; jv < vars.size(); ++jv) {
    for (size_t jj = 0; jj < nlocs_; ++jj) {
      if (mask[jj]) {
        for (size_t jl = 0; jl < nlevs_; ++jl) {
          tmp[itmp + jl] = values[ival + jl];
        }
        itmp += nlevs_;
      }
      ival += nlevs_;
    }
  }

  qg_fields_getvalsad_f90(dx.fields().toFortran(), vars, nout, locs[0], nvals, tmp[0]);
}

// -----------------------------------------------------------------------------

void InterpolatorQG::print(std::ostream & os) const {
  os << "InterpolatorQG";
}

// -----------------------------------------------------------------------------

}  // namespace qg

