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
#include <random>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"
#include "eckit/geometry/Sphere.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

#include "model/ObsVecQG.h"

namespace qg {
// -----------------------------------------------------------------------------
std::map < std::string, std::shared_ptr<ObsHelpQG> > ObsSpaceQG::theObsFileRegister_;
int ObsSpaceQG::theObsFileCount_ = 0;
// -----------------------------------------------------------------------------

ObsSpaceQG::ObsSpaceQG(const eckit::Configuration & config,
                       const util::DateTime & bgn, const util::DateTime & end)
  : oops::ObsSpaceBase(config, bgn, end), winbgn_(bgn), winend_(end), obsvars_(),
    isLocal_(false)
{
  typedef std::map< std::string, std::shared_ptr<ObsHelpQG> >::iterator otiter;

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
    helper_ = std::shared_ptr<ObsHelpQG>(new ObsHelpQG(config));
    theObsFileRegister_[ref] = helper_;
  } else {
    // File already open
    oops::Log::trace() << "ObsSpaceQG::getHelper: " << ref << " already opened." << std::endl;
    helper_ = it->second;
  }
  theObsFileCount_++;

  obsname_ = config.getString("ObsType");

  // Very UGLY!!!
  nout_ = 0;
  if (obsname_ == "Stream") nout_ = 1;
  if (obsname_ == "WSpeed") nout_ = 1;
  if (obsname_ == "Wind") nout_ = 2;
  ASSERT(nout_ > 0);
  unsigned int nvin = 0;
  if (obsname_ == "Stream") nvin = 1;
  if (obsname_ == "WSpeed") nvin = 2;
  if (obsname_ == "Wind") nvin = 2;
  ASSERT(nvin > 0);

  //  Generate locations etc... if required
  if (config.has("Generate")) {
    const eckit::LocalConfiguration gconf(config, "Generate");
    helper_->generateDistribution(gconf, obsname_, winbgn_, winend_);
  }
}

// -----------------------------------------------------------------------------

ObsSpaceQG::ObsSpaceQG(const ObsSpaceQG & obsdb,
                       const eckit::geometry::Point2 & refPoint,
                       const double & dist,
                       const int & nobs)
  : oops::ObsSpaceBase(obsdb.helper_->getConfig(), obsdb.windowStart(), obsdb.windowEnd()),
    helper_(obsdb.helper_), obsname_(obsdb.obsname_), nout_(obsdb.nout_),
    winbgn_(obsdb.winbgn_), winend_(obsdb.winend_), obsvars_(obsdb.obsvars_),
    localobs_(), isLocal_(true)
{
  oops::Log::trace() << "ObsSpaceQG for LocalObs starting" << std::endl;

  F90locs globalLocs = helper_->locations(obsdb.obsname(), winbgn_, winend_);
  for (size_t jj = 0; jj < helper_->nobs(obsname_); ++jj) {
    double lat, lon, z;
    qg_locs_element_f90(globalLocs, jj, lon, lat, z);

    eckit::geometry::Point2 obsPoint(lon, lat);
    double localDist = eckit::geometry::Sphere::distance(6.371e6, refPoint, obsPoint);
    if (localDist < dist) localobs_.push_back(jj);
  }
  oops::Log::trace() << "ObsSpaceQG for LocalObs done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsSpaceQG::getdb(const std::string & col, int & keyData) const {
  if ( isLocal_ ) {
    helper_->getdb(obsname_, col, localobs_, keyData);
  } else {
    helper_->getdb(obsname_, col, keyData);
  }
}

// -----------------------------------------------------------------------------

void ObsSpaceQG::putdb(const std::string & col, const int & keyData) const {
  ASSERT(isLocal_ == false);
  helper_->putdb(obsname_, col, keyData);
}

// -----------------------------------------------------------------------------

void ObsSpaceQG::random(const int & nn, double * xx) const {
  static util::NormalDistribution<double> dist(nn, 0.0, 1.0, getSeed());
  for (int jj = 0; jj < nn; ++jj) xx[jj] = dist[jj];
}

// -----------------------------------------------------------------------------

void ObsSpaceQG::printJo(const ObsVecQG & dy, const ObsVecQG & grad) {
  oops::Log::info() << "ObsSpaceQG::printJo not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

ObsSpaceQG::~ObsSpaceQG() {
  if ( !isLocal_ ) {
    ASSERT(theObsFileCount_ > 0);
    theObsFileCount_--;
    if (theObsFileCount_ == 0) theObsFileRegister_.clear();
  }
}

// -----------------------------------------------------------------------------
int ObsSpaceQG::nobs() const {
  if ( isLocal_ ) {
    return localobs_.size();
  } else {
    return helper_->nobs(obsname_);
  }
}
// -----------------------------------------------------------------------------

void ObsSpaceQG::print(std::ostream & os) const {
  os << "ObsSpaceQG::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace qg
