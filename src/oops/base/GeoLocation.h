/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_GEOLOCATION_H_
#define OOPS_BASE_GEOLOCATION_H_

#include <string>

#include "eckit/config/Configuration.h"
#include "oops/util/Printable.h"

namespace oops {

class GeoLocation: public util::Printable {
 public:
  GeoLocation(const double lon, const double lat):lon_(lon), lat_(lat) {}
  explicit GeoLocation(const eckit::Configuration & conf) { lon_ = conf.getDouble("lon");
                                                            lat_ = conf.getDouble("lat"); }
  ~GeoLocation() {}

  const double lon() const {return lon_;}
  const double lat() const {return lat_;}
  void getLoc(double& lon, double& lat) const {lon = lon_; lat = lat_;}

 private:
  void print(std::ostream & os) const { os << "location: " << lon_ << ", " << lat_ << std::endl; }
  double lon_, lat_;
};

}  // namespace oops

#endif  // OOPS_BASE_GEOLOCATION_H_
