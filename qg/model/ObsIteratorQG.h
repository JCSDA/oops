/*
 * (C) Copyright 2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef QG_MODEL_OBSITERATORQG_H_
#define QG_MODEL_OBSITERATORQG_H_

#include <iterator>
#include <memory>
#include <string>

#include "eckit/geometry/Point3.h"

#include "model/LocationsQG.h"

#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

namespace qg {

/// Iterator over all observations
class ObsIteratorQG:
                   public util::Printable,
                   private util::ObjectCounter<ObsIteratorQG> {
 public:
  typedef ptrdiff_t difference_type;
  typedef std::forward_iterator_tag iterator_catergory;
  typedef eckit::geometry::Point3 value_type;
  typedef eckit::geometry::Point3& reference;
  typedef eckit::geometry::Point3* pointer;

  static const std::string classname() {return "qg::ObsIteratorQG";}

  ObsIteratorQG(const ObsIteratorQG &);
  ObsIteratorQG(const LocationsQG &, int);

  bool operator==(const ObsIteratorQG &) const;
  bool operator!=(const ObsIteratorQG &) const;
  /// return location of current observation
  eckit::geometry::Point3 operator*() const;
  // pre-increment operator
  ObsIteratorQG& operator++();

 private:
  void print(std::ostream & os) const override {os << index_;}

  /// index of a current observation
  int index_;
  /// atlas field of the lons and lats
  atlas::Field locslonlat_;
};

}  // namespace qg

#endif  // QG_MODEL_OBSITERATORQG_H_
