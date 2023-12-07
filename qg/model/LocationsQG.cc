/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "atlas/array.h"
#include "atlas/field.h"
#include "atlas/functionspace/PointCloud.h"
#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "model/LocationsQG.h"
#include "model/QgFortran.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Random.h"

using atlas::array::make_datatype;
using atlas::array::make_shape;
using atlas::array::make_view;
using atlas::option::name;

namespace qg {

// -------------------------------------------------------------------------
/*! This constructor creates test locations based on settings in the config file */
LocationsQG::LocationsQG(const eckit::Configuration & config, const eckit::mpi::Comm &) {
  std::vector<double> z;
  unsigned int nrandom;

  if (config.has("lats") || config.has("nrandom")) {
    /*! These are optional lists of locations that the user can
    * specify in the config file for testing
    */
    lons_ = config.getDoubleVector("lons");
    lats_ = config.getDoubleVector("lats");
    z  = config.getDoubleVector("z");

    ASSERT(lons_.size() == lats_.size());
    ASSERT(lons_.size() == z.size());

    /*! Instead of or in addition to the specific locations in the config file
    * let the user specify a number of random locations
    */
    nrandom = config.getInt("nrandom", 0);

  } else {
    /*! If the config file does not specify otherwise, default to 4 random locations */
    nrandom = 4;
  }

  if (nrandom > 0) {
    /*! To create random locations, we need to know the dimensions of the
    * computational domain
    */
    double lonmin, lonmax, latmin, latmax, zmax;
    qg_geom_dimensions_f90(lonmin, lonmax, latmin, latmax, zmax);

    /*! Now create random locations
    * use a specified random seed for reproducibility
    */

    std::uint32_t rseed = 11;
    util::UniformDistribution<double> randlons(nrandom, lonmin, lonmax, rseed);
    util::UniformDistribution<double> randlats(nrandom, latmin, latmax);
    util::UniformDistribution<double> randzs(nrandom, 0.0, zmax);
    randlons.sort();

    for (std::size_t jj=0; jj < nrandom; ++jj) {
      lons_.push_back(randlons[jj]);
      lats_.push_back(randlats[jj]);
      z.push_back(randzs[jj]);
    }
  }

  const unsigned int nlocs = lons_.size();
  ASSERT(nlocs > 0);

  /*! render lat, lon as an atlas functionspace */
  atlas::Field field_lonlat("lonlat", make_datatype<double>(), make_shape(nlocs, 2));
  auto lonlat = make_view<double, 2>(field_lonlat);
  for ( unsigned int j = 0; j < nlocs; ++j ) {
      lonlat(j, 0) = lons_[j];
      lonlat(j, 1) = lats_[j];
  }
  pointcloud_.reset(new atlas::functionspace::PointCloud(field_lonlat));

  /*! render levels as an atlas field */
  altitude_.reset(new atlas::Field("altitude", make_datatype<double>(),
                   make_shape(nlocs)));
  auto zpc = make_view<double, 1>(*altitude_);
  for (unsigned int j = 0; j < nlocs; ++j) {
      zpc(j) = z[j];
  }

  /*! Now get times */
  if (config.has("time window")) {
    const util::DateTime winbgn(config.getString("time window.begin"));
    const util::DateTime winend(config.getString("time window.end"));
    const util::Duration window_len(winend - winbgn);
    const util::Duration dt = window_len / nlocs;

    times_.clear();
    for (unsigned int j=0; j < nlocs; ++j) {
      times_.push_back(winbgn + (j+1)*dt);
    }
  }
}

// -------------------------------------------------------------------------
/*! Deep copy constructor */
LocationsQG::LocationsQG(const LocationsQG & other) {
  pointcloud_.reset(new atlas::functionspace::PointCloud(other.lonlat()));
  altitude_.reset(new atlas::Field(*other.altitude_));
  times_ = other.times_;
  lats_ = other.lats_;
  lons_ = other.lons_;
}
// -------------------------------------------------------------------------
/*! Constructor from fields and times.  These may be obtained from the obsdb,
 * passed to C++ from Fortran
 */
LocationsQG::LocationsQG(atlas::FieldSet & fields,
                         std::vector<util::DateTime> && times) {
  pointcloud_.reset(new atlas::functionspace::PointCloud(fields.field("lonlat")));
  if (fields.has("altitude")) {
    altitude_.reset(new atlas::Field(fields.field("altitude")));
  } else {
    altitude_.reset(new atlas::Field(atlas::Field()));
  }
  times_ = times;
  lats_.resize(pointcloud_->size());
  lons_.resize(pointcloud_->size());
  atlas::Field field_lonlat = pointcloud_->lonlat();
  auto lonlat = make_view<double, 2>(field_lonlat);
  for (atlas::idx_t jloc = 0; jloc < pointcloud_->size(); ++jloc) {
    lats_[jloc] = lonlat(jloc, 1);
    lons_[jloc] = lonlat(jloc, 0);
  }
}
// -------------------------------------------------------------------------
void LocationsQG::print(std::ostream & os) const {
  int nobs = pointcloud_->size();
  atlas::Field field_lonlat = pointcloud_->lonlat();
  auto lonlat = make_view<double, 2>(field_lonlat);
  auto z = make_view<double, 1>(*altitude_);
  std::ios_base::fmtflags fmt = os.flags();
  std::streamsize ss = os.precision();
  for (size_t jj=0; jj < static_cast<size_t>(nobs); ++jj) {
    os << "location " << jj << std::fixed << std::setprecision(2) << ": lon = " << lonlat(jj, 0)
       << ", lat = " << lonlat(jj, 1) << ", z = " << z(jj) << std::endl;
  }
  os.setf(fmt);
  os.precision(ss);
}
// -------------------------------------------------------------------------
}  // namespace qg
