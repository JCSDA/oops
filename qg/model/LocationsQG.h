/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_LOCATIONSQG_H_
#define QG_MODEL_LOCATIONSQG_H_

#include <iomanip>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/field.h"
#include "atlas/functionspace/PointCloud.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/mpi/Comm.h"

#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"

#include "oops/qg/QgFortran.h"

namespace qg {

// -----------------------------------------------------------------------------
/// implements the oops SampledLocations interface for the QG obs type.
class LocationsQG : public util::Printable,
                    private util::ObjectCounter<LocationsQG> {
 public:
  static const std::string classname() {return "qg::LocationsQG";}

// Constructors and basic operators
  LocationsQG(atlas::FieldSet &, std::vector<util::DateTime> &&);
  LocationsQG(const eckit::Configuration &, const eckit::mpi::Comm &);
  LocationsQG(const LocationsQG &);
  ~LocationsQG() {}

// Utilities
  int size() const {return pointcloud_->size();}
  atlas::functionspace::PointCloud & pointcloud() {return *pointcloud_;}
  atlas::Field lonlat() const {return pointcloud_->lonlat();}
  atlas::Field & altitude() {ASSERT(altitude_); return *altitude_;}
  const std::vector<double> & latitudes() const {return lats_;}
  const std::vector<double> & longitudes() const {return lons_;}
  const std::vector<util::DateTime> & times() const {return times_;}
  util::DateTime & times(size_t idx) {return times_[idx];}

 private:
  void print(std::ostream &) const;
  std::unique_ptr<atlas::functionspace::PointCloud> pointcloud_;
  std::unique_ptr<atlas::Field> altitude_;
  std::vector<util::DateTime> times_;
  std::vector<double> lats_;
  std::vector<double> lons_;
};

}  // namespace qg

#endif  // QG_MODEL_LOCATIONSQG_H_
