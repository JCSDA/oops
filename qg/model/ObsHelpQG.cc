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

#include "model/ObsHelpQG.h"

#include <string>

#include "eckit/config/Configuration.h"
#include "model/QgFortran.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace qg {

// -----------------------------------------------------------------------------

ObsHelpQG::ObsHelpQG(const eckit::Configuration & config)
  : config_(config) {
  const eckit::Configuration * configc = &config;
  qg_obsdb_setup_f90(keyHelp_, &configc);
  oops::Log::trace() << "ObsHelpQG constructed" << std::endl;
}

// -----------------------------------------------------------------------------

ObsHelpQG::~ObsHelpQG() {
  qg_obsdb_delete_f90(keyHelp_);
  oops::Log::trace() << "ObsHelpQG destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsHelpQG::putdb(const std::string & obsname, const std::string & col, const int & keyFvec) {
  oops::Log::trace() << "ObsHelpQG:putdb obsname = " << obsname << ", col = " << col << std::endl;
  qg_obsdb_put_f90(keyHelp_, obsname.size(), obsname.c_str(), col.size(), col.c_str(), keyFvec);
}

// -----------------------------------------------------------------------------

void ObsHelpQG::getdb(const std::string & obsname, const std::string & col, int & keyFvec) const {
  oops::Log::trace() << "ObsHelpQG:getdb obsname = " << obsname << ", col = " << col << std::endl;
  qg_obsdb_get_f90(keyHelp_, obsname.size(), obsname.c_str(), col.size(), col.c_str(), keyFvec);
  oops::Log::trace() << "ObsHelpQG:getdb done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsHelpQG::getdb(const std::string & obsname, const std::string & col,
                      const std::vector<int> & localidx, int & keyFvec) const {
  oops::Log::trace() << "ObsHelpQG:getdb for LocalObs obsname = " << obsname
                     << ", col = " << col << std::endl;
  qg_obsdb_get_local_f90(keyHelp_, obsname.size(), obsname.c_str(), col.size(),
                         col.c_str(), localidx.size(), localidx.data(), keyFvec);
  oops::Log::trace() << "ObsHelpQG:getdb for LocalObs done" << std::endl;
}

// -----------------------------------------------------------------------------

bool ObsHelpQG::has(const std::string & obsname, const std::string & col) const {
  int ii;
  qg_obsdb_has_f90(keyHelp_, obsname.size(), obsname.c_str(), col.size(), col.c_str(), ii);
  return ii;
}

// -----------------------------------------------------------------------------

F90locs ObsHelpQG::locations(const std::string & obsname,
                             const util::DateTime & t1, const util::DateTime & t2) const {
  const util::DateTime * p1 = &t1;
  const util::DateTime * p2 = &t2;
  F90locs key_locs;
  qg_obsdb_locations_f90(keyHelp_, obsname.size(), obsname.c_str(), &p1, &p2, key_locs);
  return key_locs;
}

// -----------------------------------------------------------------------------

void ObsHelpQG::generateDistribution(const eckit::Configuration & config,
                                     const std::string & obsname,
                                     const util::DateTime & t1,
                                     const util::DateTime & t2) {
  const eckit::Configuration * configc = &config;
  const util::Duration first(config.getString("begin"));
  const util::DateTime start(t1 + first);
  const util::Duration freq(config.getString("obs_period"));
  int nobstimes = 0;
  util::DateTime now(start);
  while (now <= t2) {
    ++nobstimes;
    now += freq;
  }
  const util::DateTime * bgn = &start;
  const util::Duration * stp = &freq;
  int iobs;
  qg_obsdb_generate_f90(keyHelp_, obsname.size(), obsname.c_str(), &configc,
                        &bgn, &stp, nobstimes, iobs);
}

// -----------------------------------------------------------------------------

unsigned int ObsHelpQG::nobs(const std::string & obsname) const {
  int iobs;
  qg_obsdb_nobs_f90(keyHelp_, obsname.size(), obsname.c_str(), iobs);
  return iobs;
}

// -----------------------------------------------------------------------------


}  // namespace qg
