/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include <utility>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"

#include "model/ObsIteratorQG.h"

// -----------------------------------------------------------------------------
namespace qg {

ObsIteratorQG::ObsIteratorQG(const ObsIteratorQG & other): index_(other.index_),
      locslonlat_(other.locslonlat_) {}

// -----------------------------------------------------------------------------
ObsIteratorQG::ObsIteratorQG(const LocationsQG & locations, int index):
      index_(index), locslonlat_(locations.lonlat()) {}

// -----------------------------------------------------------------------------
bool ObsIteratorQG::operator==(const ObsIteratorQG & other) const {
  return (index_ == other.index_);
}

// -----------------------------------------------------------------------------
bool ObsIteratorQG::operator!=(const ObsIteratorQG & other) const {
  return (index_!= other.index_);
}

// -----------------------------------------------------------------------------
eckit::geometry::Point2 ObsIteratorQG::operator*() const {
  auto lonlat = atlas::array::make_view<double, 2>(locslonlat_);
  return eckit::geometry::Point2(lonlat(index_, 0), lonlat(index_, 1));
}

// -----------------------------------------------------------------------------
ObsIteratorQG& ObsIteratorQG::operator++() {
  index_++;
  return *this;
}

// -----------------------------------------------------------------------------

}  // namespace qg
