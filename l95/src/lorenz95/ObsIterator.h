/*
 * (C) Copyright 2021 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef LORENZ95_OBSITERATOR_H_
#define LORENZ95_OBSITERATOR_H_

#include <iterator>
#include <string>
#include <vector>

#include "eckit/geometry/Point3.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace lorenz95 {

/// Iterator over all observations
class ObsIterator: public util::Printable,
                   private util::ObjectCounter<ObsIterator> {
 public:
  // for std::iterator_traits
  typedef std::forward_iterator_tag iterator_category;
  typedef eckit::geometry::Point3 value_type;
  typedef eckit::geometry::Point3& reference;
  typedef eckit::geometry::Point3* pointer;
  typedef std::ptrdiff_t difference_type;

  static const std::string classname() {return "lorenz95::ObsIterator";}

  ObsIterator(const std::vector<double> & locations, int index);

  bool operator==(const ObsIterator &) const;
  bool operator!=(const ObsIterator &) const;
  /// return location of current observation
  eckit::geometry::Point3 operator*() const;
  // pre-increment operator
  ObsIterator& operator++();

 private:
  void print(std::ostream & os) const override {os << index_;}

  /// locations (in 1D) of all observations
  const std::vector<double> locations_;
  /// index of a current observation
  int index_;
};

}  // namespace lorenz95

#endif  // LORENZ95_OBSITERATOR_H_
