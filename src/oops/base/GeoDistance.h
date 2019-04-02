/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_GEODISTANCE_H_
#define OOPS_BASE_GEODISTANCE_H_

#include <string>

#include "eckit/config/Configuration.h"
#include "oops/util/Printable.h"

namespace oops {

class GeoDistance: public util::Printable {
 public:
  explicit GeoDistance(const double & dist): dist_(dist) {}
  explicit GeoDistance(const eckit::Configuration & conf) { dist_ = conf.getDouble("distance"); }
  ~GeoDistance() {}

  double distance() const {return dist_;}

 private:
  void print(std::ostream & os) const { os << "distance: " << dist_ << std::endl; }
  double dist_;
};

}  // namespace oops

#endif  // OOPS_BASE_GEODISTANCE_H_
