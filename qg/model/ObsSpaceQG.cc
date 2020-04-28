/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/ObsSpaceQG.h"

#include <map>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/geometry/Sphere.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "model/ObsVecQG.h"

namespace qg {
// -----------------------------------------------------------------------------
// initialization for the static map
std::map < std::string, F90odb > ObsSpaceQG::theObsFileRegister_;
int ObsSpaceQG::theObsFileCount_ = 0;

// -----------------------------------------------------------------------------

ObsSpaceQG::ObsSpaceQG(const eckit::Configuration & config, const eckit::mpi::Comm & comm,
                       const util::DateTime & bgn, const util::DateTime & end)
  : oops::ObsSpaceBase(config, comm, bgn, end), obsname_(config.getString("ObsType")),
    winbgn_(bgn), winend_(end), obsvars_(), isLocal_(false), comm_(comm)
{
  typedef std::map< std::string, F90odb >::iterator otiter;

  std::string ofin("-");
  if (config.has("ObsDataIn")) {
    ofin = config.getString("ObsDataIn.obsfile");
  }
  std::string ofout("-");
  if (config.has("ObsDataOut")) {
    ofout = config.getString("ObsDataOut.obsfile");
  }
  oops::Log::trace() << "ObsSpaceQG: Obs files are: " << ofin << " and " << ofout << std::endl;
  std::string ref = ofin + ofout;
  if (ref == "--") {
    ABORT("Underspecified observation files.");
  }

  otiter it = theObsFileRegister_.find(ref);
  if ( it == theObsFileRegister_.end() ) {
    // Open new file
    oops::Log::trace() << "ObsSpaceQG::getHelper: " << "Opening " << ref << std::endl;
    qg_obsdb_setup_f90(key_, config);
    theObsFileRegister_[ref] = key_;
  } else {
    // File already open
    oops::Log::trace() << "ObsSpaceQG::getHelper: " << ref << " already opened." << std::endl;
    key_ = it->second;
  }
  theObsFileCount_++;

  // Set variables simulated for different obstypes
  if (obsname_ == "Stream") obsvars_.push_back("Stream");
  if (obsname_ == "WSpeed") obsvars_.push_back("WSpeed");
  if (obsname_ == "Wind") {
    obsvars_.push_back("Uwind");
    obsvars_.push_back("Vwind");
  }

  //  Generate locations etc... if required
  if (config.has("Generate")) {
    const eckit::LocalConfiguration gconf(config, "Generate");
    const util::Duration first(gconf.getString("begin"));
    const util::DateTime start(winbgn_ + first);
    const util::Duration freq(gconf.getString("obs_period"));
    int nobstimes = 0;
    util::DateTime now(start);
    while (now <= winend_) {
      ++nobstimes;
      now += freq;
    }
    int iobs;
    qg_obsdb_generate_f90(key_, obsname_.size(), obsname_.c_str(), gconf,
                          start, freq, nobstimes, iobs);
  }
}

// -----------------------------------------------------------------------------

ObsSpaceQG::ObsSpaceQG(const ObsSpaceQG & obsdb,
                       const eckit::geometry::Point2 & refPoint,
                       const eckit::Configuration & conf)
  : oops::ObsSpaceBase(eckit::LocalConfiguration(), obsdb.comm_,
                       obsdb.windowStart(), obsdb.windowEnd()),
    key_(obsdb.key_), obsname_(obsdb.obsname_),
    winbgn_(obsdb.winbgn_), winend_(obsdb.winend_), obsvars_(obsdb.obsvars_),
    localobs_(), isLocal_(true), comm_(obsdb.comm_)
{
  oops::Log::trace() << "ObsSpaceQG for LocalObs starting" << std::endl;
  const double dist = conf.getDouble("lengthscale");
  F90locs key_locs;
  int iobs;
  // get number of obs
  qg_obsdb_nobs_f90(key_, obsname_.size(), obsname_.c_str(), iobs);
  // get locations of all obs
  qg_obsdb_locations_f90(key_, obsname_.size(), obsname_.c_str(), winbgn_, winend_, key_locs);
  for (size_t jj = 0; jj < iobs; ++jj) {
    double lat, lon, z;
    qg_locs_element_f90(key_locs, jj, lon, lat, z);
    eckit::geometry::Point2 obsPoint(lon, lat);
    double localDist = eckit::geometry::Sphere::distance(6.371e6, refPoint, obsPoint);
    if (localDist < dist) localobs_.push_back(jj);
  }
  oops::Log::trace() << "ObsSpaceQG for LocalObs done" << std::endl;
}

// -----------------------------------------------------------------------------

ObsSpaceQG::~ObsSpaceQG() {
  if ( !isLocal_ ) {
    ASSERT(theObsFileCount_ > 0);
    theObsFileCount_--;
    if (theObsFileCount_ == 0) {
      theObsFileRegister_.clear();
      qg_obsdb_delete_f90(key_);
    }
  }
}

// -----------------------------------------------------------------------------

void ObsSpaceQG::getdb(const std::string & col, int & keyData) const {
  if ( isLocal_ ) {
    qg_obsdb_get_local_f90(key_, obsname_.size(), obsname_.c_str(), col.size(),
                         col.c_str(), localobs_.size(), localobs_.data(), keyData);
  } else {
    qg_obsdb_get_f90(key_, obsname_.size(), obsname_.c_str(), col.size(), col.c_str(), keyData);
  }
}

// -----------------------------------------------------------------------------

void ObsSpaceQG::putdb(const std::string & col, const int & keyData) const {
  // not implemented for local ObsSpace
  ASSERT(isLocal_ == false);
  qg_obsdb_put_f90(key_, obsname_.size(), obsname_.c_str(), col.size(), col.c_str(), keyData);
}

// -----------------------------------------------------------------------------

bool ObsSpaceQG::has(const std::string & col) const {
  int ii;
  qg_obsdb_has_f90(key_, obsname_.size(), obsname_.c_str(), col.size(), col.c_str(), ii);
  return ii;
}

// -----------------------------------------------------------------------------

LocationsQG * ObsSpaceQG::locations(const util::DateTime & t1, const util::DateTime & t2) const {
  F90locs key_locs;
  qg_obsdb_locations_f90(key_, obsname_.size(), obsname_.c_str(), t1, t2, key_locs);
  return new LocationsQG(key_locs);
}

// -----------------------------------------------------------------------------

void ObsSpaceQG::printJo(const ObsVecQG & dy, const ObsVecQG & grad) const {
  oops::Log::info() << "ObsSpaceQG::printJo not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

int ObsSpaceQG::nobs() const {
  if ( isLocal_ ) {
    return localobs_.size();
  } else {
    int iobs;
    qg_obsdb_nobs_f90(key_, obsname_.size(), obsname_.c_str(), iobs);
    return iobs;
  }
}
// -----------------------------------------------------------------------------

void ObsSpaceQG::print(std::ostream & os) const {
  os << "ObsSpace for " << obsname_ << ", " << this->nobs() << " obs" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace qg
