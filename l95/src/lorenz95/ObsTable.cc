/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/ObsTable.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "util/abor1_cpp.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/Logger.h"

#include "lorenz95/ObsVec1D.h"

using std::string;
using std::endl;
using std::ifstream;
using std::ofstream;


using oops::Log;

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------

ObsTable::ObsTable(const eckit::Configuration & config,
                   const util::DateTime & bgn, const util::DateTime & end)
  : winbgn_(bgn), winend_(end)
{
  nameIn_.clear();
  nameOut_.clear();
  if (config.has("ObsData")) {
    const eckit::LocalConfiguration dataConfig(config, "ObsData");
    if (dataConfig.has("ObsDataIn")) {
      nameIn_ = dataConfig.getString("ObsDataIn.filename");
      Log::trace() << "ObsTable::ObsTable reading observations from " << nameIn_ << std::endl;
      otOpen(nameIn_);
    }
    if (dataConfig.has("ObsDataOut")) {
      nameOut_ = dataConfig.getString("ObsDataOut.filename");
    }
  }
  Log::trace() << "ObsTable::ObsTable created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsTable::~ObsTable() {
  if (!nameOut_.empty()) {
    Log::trace() << "ObsTable::~ObsTable saving nameOut = " << nameOut_ << std::endl;
    otWrite(nameOut_);
  }
  Log::trace() << "ObsTable::ObsTable destructed" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTable::putdb(const string & col, const std::vector<double> & vec) const {
  ASSERT(vec.size() == nobs());
  ASSERT(data_.find(col) == data_.end());
  data_.insert(std::pair<string, std::vector<double> >(col, vec));
}

// -----------------------------------------------------------------------------

void ObsTable::getdb(const string & col, std::vector<double> & vec) const {
  std::map<string, std::vector<double> >::const_iterator ic = data_.find(col);
  ASSERT(ic != data_.end());
  vec.resize(nobs());
  for (unsigned int jobs = 0; jobs < nobs(); ++jobs) {
    vec[jobs] = ic->second[jobs];
  }
}

// -----------------------------------------------------------------------------

std::vector<double> ObsTable::locations(const util::DateTime & t1,
                                        const util::DateTime & t2) const {
  std::vector<int> olist = timeSelect(t1, t2);
  const int nobs = olist.size();
  std::vector<double> locs(nobs);
  for (int jobs = 0; jobs < nobs; ++jobs) locs[jobs]=locations_[olist[jobs]];
  return locs;
}

// -----------------------------------------------------------------------------

std::vector<int> ObsTable::timeSelect(const util::DateTime & t1,
                                      const util::DateTime & t2) const {
  std::vector<int> mask;
  for (unsigned int jobs = 0; jobs < nobs(); ++jobs)
    if (times_[jobs] > t1 && times_[jobs] <= t2) mask.push_back(jobs);
  return mask;
}

// -----------------------------------------------------------------------------

void ObsTable::generateDistribution(const eckit::Configuration & config) {
  Log::trace() << "ObsTable::generateDistribution starting" << std::endl;

  util::Duration first(config.getString("begin"));
  util::Duration last(winend_-winbgn_);
  if (config.has("end")) {
    last = util::Duration(config.getString("end"));
  }
  util::Duration freq(config.getString("obs_frequency"));

  int nobstimes = 0;
  util::Duration step(first);
  while (step <= last) {
    ++nobstimes;
    step += freq;
  }

  const unsigned int nobs_locations = config.getInt("obs_density");
  const unsigned int nobs = nobs_locations*nobstimes;
  double dx = 1.0/static_cast<double>(nobs_locations);
  Log::trace() << "ObservationL95:generateDistribution nobs=" << nobs << std::endl;

  times_.resize(nobs);
  locations_.resize(nobs);

  unsigned int iobs = 0;
  util::DateTime now(winbgn_);
  step = first;
  while (step <= last) {
    now = winbgn_ + step;
    for (unsigned int jobs = 0; jobs < nobs_locations; ++jobs) {
      double xpos = jobs*dx;
      times_[iobs] = now;
      locations_[iobs] = xpos;
      ++iobs;
    }
    step += freq;
  }
  ASSERT(iobs == nobs);

  Log::trace() << "ObsTable::generateDistribution done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTable::printJo(const ObsVec1D & ydep, const ObsVec1D & grad) {
  Log::info() << "ObsTable::printJo not implemented" << std::endl;
}

// -----------------------------------------------------------------------------
//  ObsTable Private Methods
// -----------------------------------------------------------------------------

void ObsTable::otOpen(const string & filename) {
  Log::trace() << "ObsTable::ot_read starting" << std::endl;
  ifstream fin(filename.c_str());
  if (!fin.is_open()) ABORT("ObsTable::otOpen: Error opening file");

  int ncol, nobs;
  fin >> ncol;

  std::vector<string> colnames;
  for (int jc = 0; jc < ncol; ++jc) {
    string col;
    fin >> col;
    colnames.push_back(col);
  }

  fin >> nobs;
  locations_.resize(nobs);
  std::vector<double> newcol(nobs);
  for (int jc = 0; jc < ncol; ++jc) {
    ASSERT(data_.find(colnames[jc]) == data_.end());
    data_.insert(std::pair<string, std::vector<double> >(colnames[jc], newcol));
  }

  times_.clear();
  int jjj;
  for (int jobs = 0; jobs < nobs; ++jobs) {
    fin >> jjj;
    ASSERT(jjj == jobs);
    string sss;
    fin >> sss;
    util::DateTime ttt(sss);
    times_.push_back(ttt);
    fin >> locations_[jobs];
    for (std::map<string, std::vector<double> >::iterator jo = data_.begin();
         jo != data_.end(); ++jo)
      fin >> jo->second[jobs];
  }

  fin.close();
  Log::trace() << "ObsTable::ot_read done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTable::otWrite(const string & filename) const {
  Log::trace() << "ObsTable::otWrite writing " << filename << std::endl;
  ofstream fout(filename.c_str());
  if (!fout.is_open()) ABORT("ObsTable::otWrite: Error opening file");

  int ncol = data_.size();
  fout << ncol << endl;

  for (std::map<string, std::vector<double> >::const_iterator jo = data_.begin();
       jo != data_.end(); ++jo)
    fout << jo->first << endl;

  int nobs = times_.size();
  fout << nobs << endl;
  for (int jobs = 0; jobs < nobs; ++jobs) {
    fout << jobs;
    fout << "  " << times_[jobs];
    fout << "  " << locations_[jobs];
    for (std::map<string, std::vector<double> >::const_iterator jo = data_.begin();
         jo != data_.end(); ++jo)
      fout << "  " << jo->second[jobs];
    fout << endl;
  }

  fout.close();
  Log::trace() << "ObsTable::otWrite done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTable::print(std::ostream & os) const {
  os << "ObsTable::print not implemented";
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
