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
  : grid_(grid), nlocs_(lats.size()), locs_(2 * nlocs_)
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
  const size_t nlevs = grid_.levels();
  ASSERT(mask.size() == nlocs_);
  ASSERT(values.size() == vars.size() * nlevs * nlocs_);

  size_t nout = 0;
  std::vector<double> locs;
  for (size_t jj = 0; jj < nlocs_; ++jj) {
    if (mask[jj]) {
      locs.push_back(locs_[2*jj]);
      locs.push_back(locs_[2*jj+1]);
      ++nout;
    }
  }
  const size_t nvals = vars.size() * nlevs * nout;
  std::vector<double> tmp(nvals);

  qg_fields_getvals_f90(flds.toFortran(), vars, nout, locs[0], nvals, tmp[0]);

  // Reorder data: tmp is levels-varies-fastest, want locs-varies-fastest
  size_t ival = 0;
  for (size_t jv = 0; jv < vars.size(); ++jv) {
    const size_t jvStart = jv * nlevs * nout;
    for (size_t jl = 0; jl < nlevs; ++jl) {
      size_t itmp = jvStart + jl;
      for (size_t jj = 0; jj < nlocs_; ++jj) {
        if (mask[jj]) {
          values[ival] = tmp[itmp];
          itmp += nlevs;
        }
        ++ival;
      }
    }
  }
}

// -----------------------------------------------------------------------------

void InterpolatorQG::applyAD(const oops::Variables & vars, IncrementQG & dx,
                             const std::vector<bool> & mask,
                             const std::vector<double> & values) const {
  const size_t nlevs = grid_.levels();
  ASSERT(mask.size() == nlocs_);
  ASSERT(values.size() == vars.size() * nlevs * nlocs_);

  size_t nout = 0;
  std::vector<double> locs;
  for (size_t jj = 0; jj < nlocs_; ++jj) {
    if (mask[jj]) {
      locs.push_back(locs_[2*jj]);
      locs.push_back(locs_[2*jj+1]);
      ++nout;
    }
  }
  const size_t nvals = vars.size() * nlevs * nout;
  std::vector<double> tmp(nvals);

  // (Adjoint of) Reorder data: tmp is levels-varies-fastest, want locs-varies-fastest
  size_t ival = 0;
  for (size_t jv = 0; jv < vars.size(); ++jv) {
    const size_t jvStart = jv * nlevs * nout;
    for (size_t jl = 0; jl < nlevs; ++jl) {
      size_t itmp = jvStart + jl;
      for (size_t jj = 0; jj < nlocs_; ++jj) {
        if (mask[jj]) {
          tmp[itmp] = values[ival];
          itmp += nlevs;
        }
        ++ival;
      }
    }
  }

  qg_fields_getvalsad_f90(dx.fields().toFortran(), vars, nout, locs[0], nvals, tmp[0]);
}

// -----------------------------------------------------------------------------

void InterpolatorQG::apply(const oops::Variables & vars, const atlas::FieldSet & fset,
                           const std::vector<bool> & mask,
                           std::vector<double> & values) const {
  const size_t nlevs = grid_.levels();
  ASSERT(mask.size() == nlocs_);
  ASSERT(values.size() == vars.size() * nlevs * nlocs_);

  size_t nout = 0;
  std::vector<double> locs;
  for (size_t jj = 0; jj < nlocs_; ++jj) {
    if (mask[jj]) {
      locs.push_back(locs_[2*jj]);
      locs.push_back(locs_[2*jj+1]);
      ++nout;
    }
  }
  const size_t nvals = vars.size() * nlevs * nout;
  std::vector<double> tmp(nvals);

  qg_getvalues_interp_f90(grid_.toFortran(), fset.get(), vars, nout, locs[0], nvals, tmp[0]);

  // Reorder data: tmp is levels-varies-fastest, want locs-varies-fastest
  size_t ival = 0;
  for (size_t jv = 0; jv < vars.size(); ++jv) {
    const size_t jvStart = jv * nlevs * nout;
    for (size_t jl = 0; jl < nlevs; ++jl) {
      size_t itmp = jvStart + jl;
      for (size_t jj = 0; jj < nlocs_; ++jj) {
        if (mask[jj]) {
          values[ival] = tmp[itmp];
          itmp += nlevs;
        }
        ++ival;
      }
    }
  }
}

// -----------------------------------------------------------------------------

void InterpolatorQG::applyAD(const oops::Variables & vars, atlas::FieldSet & fset,
                             const std::vector<bool> & mask,
                             const std::vector<double> & values) const {
  const size_t nlevs = grid_.levels();
  ASSERT(mask.size() == nlocs_);
  ASSERT(values.size() == vars.size() * nlevs * nlocs_);

  size_t nout = 0;
  std::vector<double> locs;
  for (size_t jj = 0; jj < nlocs_; ++jj) {
    if (mask[jj]) {
      locs.push_back(locs_[2*jj]);
      locs.push_back(locs_[2*jj+1]);
      ++nout;
    }
  }
  const size_t nvals = vars.size() * nlevs * nout;
  std::vector<double> tmp(nvals);

  // (Adjoint of) Reorder data: tmp is levels-varies-fastest, want locs-varies-fastest
  size_t ival = 0;
  for (size_t jv = 0; jv < vars.size(); ++jv) {
    const size_t jvStart = jv * nlevs * nout;
    for (size_t jl = 0; jl < nlevs; ++jl) {
      size_t itmp = jvStart + jl;
      for (size_t jj = 0; jj < nlocs_; ++jj) {
        if (mask[jj]) {
          tmp[itmp] = values[ival];
          itmp += nlevs;
        }
        ++ival;
      }
    }
  }

  qg_getvalues_interp_ad_f90(grid_.toFortran(), fset.get(), vars, nout, locs[0], nvals, tmp[0]);
}

// -----------------------------------------------------------------------------

void InterpolatorQG::print(std::ostream & os) const {
  os << "InterpolatorQG";
}

// -----------------------------------------------------------------------------

}  // namespace qg

