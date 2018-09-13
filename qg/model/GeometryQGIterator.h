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

#include "model/GeometryQG.h"
#include "model/QgFortran.h"

#include "oops/base/GeoLocation.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace oops {
  class GeoLocation;
}

namespace qg {

class GeometryQG;

// -----------------------------------------------------------------------------
class GeometryQGIterator: public std::iterator<std::forward_iterator_tag,
                                               oops::GeoLocation>,
                          public util::Printable,
                          private util::ObjectCounter<GeometryQGIterator> {
 public:
  static const std::string classname() {return "qg::GeometryQGIterator";}

  GeometryQGIterator(const GeometryQGIterator &);
  explicit GeometryQGIterator(const GeometryQG & geom, const int & index = 1);
  ~GeometryQGIterator();

  bool operator==(const GeometryQGIterator &) const;
  bool operator!=(const GeometryQGIterator &) const;
  oops::GeoLocation operator*() const;
  GeometryQGIterator& operator++();

  F90iter & toFortran() {return keyIter_;}
  const F90iter & toFortran() const {return keyIter_;}

 private:
  void print(std::ostream &) const;
  F90iter keyIter_;
};

}  // namespace qg

#endif  // QG_MODEL_GEOMETRYQGITERATOR_H_
