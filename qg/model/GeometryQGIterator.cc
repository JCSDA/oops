/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "eckit/config/Configuration.h"
#include "model/GeometryQGIterator.h"
#include "model/QgFortran.h"
#include "oops/util/Logger.h"

// -----------------------------------------------------------------------------

namespace qg {

// -----------------------------------------------------------------------------

GeometryQGIterator::GeometryQGIterator(const GeometryQGIterator& iter) {
  qg_geom_iter_clone_f90(keyIter_, iter.toFortran());
}

// -----------------------------------------------------------------------------

GeometryQGIterator::GeometryQGIterator(const GeometryQG& geom, const int & index) {
  qg_geom_iter_setup_f90(keyIter_, geom.toFortran(), index);
}


// -----------------------------------------------------------------------------

GeometryQGIterator::~GeometryQGIterator() {
  qg_geom_iter_delete_f90(keyIter_);
}

// -----------------------------------------------------------------------------

bool GeometryQGIterator::operator==(const GeometryQGIterator & other) const {
  int equals = 0;
  qg_geom_iter_equals_f90(keyIter_, other.toFortran(), equals);
  return (equals == 1);
}

// -----------------------------------------------------------------------------

bool GeometryQGIterator::operator!=(const GeometryQGIterator & other) const {
  int equals = 0;
  qg_geom_iter_equals_f90(keyIter_, other.toFortran(), equals);
  return (equals == 0);
}

// -----------------------------------------------------------------------------

eckit::geometry::Point3 GeometryQGIterator::operator*() const {
  double lat, lon;
  qg_geom_iter_current_f90(keyIter_, lat, lon);
  return eckit::geometry::Point3(lat, lon, 0.0);
}

// -----------------------------------------------------------------------------

GeometryQGIterator& GeometryQGIterator::operator++() {
  qg_geom_iter_next_f90(keyIter_);
  return *this;
}

// -----------------------------------------------------------------------------

void GeometryQGIterator::print(std::ostream & os) const {
  os << "GeometryQGIterator, key: " <<  keyIter_;
}

// -----------------------------------------------------------------------------

}  // namespace qg
