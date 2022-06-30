/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_ITERATOR_H_
#define LORENZ95_ITERATOR_H_

#include <iterator>
#include <string>
#include <vector>

#include "eckit/geometry/Point3.h"

#include "lorenz95/Resolution.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace lorenz95 {

class Resolution;

// -----------------------------------------------------------------------------
class Iterator: public std::iterator<std::forward_iterator_tag,
                                     eckit::geometry::Point3>,
                public util::Printable,
                private util::ObjectCounter<Iterator> {
 public:
  static const std::string classname() {return "lorenz95::Iterator";}

  Iterator(const Resolution & res, const int & index);

  bool operator==(const Iterator &) const;
  bool operator!=(const Iterator &) const;
  eckit::geometry::Point3 operator*() const;
  Iterator& operator++();

  int index() const {return index_;}

 private:
  void print(std::ostream & os) const override {os << index_;}
  const int res_;
  int index_;
};

}  // namespace lorenz95

#endif  // LORENZ95_ITERATOR_H_
