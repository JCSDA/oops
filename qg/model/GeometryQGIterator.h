/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_GEOMETRYQGITERATOR_H_
#define QG_MODEL_GEOMETRYQGITERATOR_H_

#include <iterator>
#include <string>

#include "eckit/geometry/Point3.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "oops/qg/GeometryQG.h"
#include "oops/qg/QgFortran.h"

namespace qg {

class GeometryQG;

// -----------------------------------------------------------------------------
class GeometryQGIterator: public util::Printable,
                          private util::ObjectCounter<GeometryQGIterator> {
 public:
  typedef std::forward_iterator_tag iterator_category;
  typedef eckit::geometry::Point3 value_type;
  typedef eckit::geometry::Point3& reference;
  typedef eckit::geometry::Point3* pointer;
  typedef ptrdiff_t difference_type;

  static const std::string classname() {return "qg::GeometryQGIterator";}

  GeometryQGIterator(const GeometryQGIterator &);
  explicit GeometryQGIterator(const GeometryQG & geom, const int & index = 1);
  ~GeometryQGIterator();

  bool operator==(const GeometryQGIterator &) const;
  bool operator!=(const GeometryQGIterator &) const;
  eckit::geometry::Point3 operator*() const;
  // pre-increment operator
  GeometryQGIterator& operator++();

  const F90iter & toFortran() const {return keyIter_;}

 private:
  void print(std::ostream &) const;
  F90iter keyIter_;
};

}  // namespace qg

#endif  // QG_MODEL_GEOMETRYQGITERATOR_H_
