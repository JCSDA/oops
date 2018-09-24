/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_GEOLOCATION_H_
#define OOPS_BASE_GEOLOCATION_H_

#include <string>

#include "oops/util/Printable.h"

namespace oops {

class GeoLocation: public util::Printable {
 public:
  GeoLocation(const double lon, const double lat):lon_(lon), lat_(lat) {}
  ~GeoLocation() {}

  void getLoc(double& lon, double& lat) const {lon = lon_; lat = lat_;}

 private:
  void print(std::ostream & os) const { os << "location: " << lon_ << ", " << lat_ << std::endl; }
  double lon_, lat_;
};

}  // namespace oops

#endif  // OOPS_BASE_GEOLOCATION_H_
