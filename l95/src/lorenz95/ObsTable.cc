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
#include <cmath>
#include <fstream>
#include <iostream>
#include <map>
#include <memory>
#include <string>
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "lorenz95/ObsIterator.h"
#include "oops/mpi/mpi.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"
#include "oops/util/Random.h"
#include "oops/util/stringFunctions.h"

namespace sf = util::stringfunctions;

// -----------------------------------------------------------------------------
namespace lorenz95 {
// -----------------------------------------------------------------------------

ObsTable::ObsTable(const eckit::Configuration & config, const eckit::mpi::Comm & comm,
                   const util::DateTime & bgn, const util::DateTime & end,
                   const eckit::mpi::Comm & timeComm)
  : oops::ObsSpaceBase(config, comm, bgn, end), winbgn_(bgn), winend_(end), comm_(timeComm),
    obsvars_()
{
  oops::Log::trace() << "ObsTable::ObsTable starting" << std::endl;
  nameIn_.clear();
  nameOut_.clear();
  if (config.has("obsdatain")) {
    nameIn_ = config.getString("obsdatain");
    otOpen(nameIn_);
  }
  //  Generate locations etc... if required
  if (config.has("generate")) {
    const eckit::LocalConfiguration gconf(config, "generate");
    generateDistribution(gconf);
  }
  if (config.has("obsdataout")) {
    nameOut_ = config.getString("obsdataout");
    sf::swapNameMember(config, nameOut_);
  }
  oops::Log::trace() << "ObsTable::ObsTable created nobs = " << nobs() << std::endl;
}

// -----------------------------------------------------------------------------

ObsTable::~ObsTable() {
  if (!nameOut_.empty()) {
    otWrite(nameOut_);
  }
  oops::Log::trace() << "ObsTable::ObsTable destructed" << std::endl;
}

// -----------------------------------------------------------------------------

bool ObsTable::has(const std::string & col) const {
  return (data_.find(col) != data_.end());
}

// -----------------------------------------------------------------------------

void ObsTable::putdb(const std::string & col, const std::vector<int> & vec) const {
  std::vector<double> tmp(vec.size());
  int intmiss = util::missingValue(intmiss);
  double doublemiss = util::missingValue(doublemiss);
  for (size_t jobs = 0; jobs < vec.size(); ++jobs) {
    if (vec[jobs] == intmiss) {
      tmp[jobs] = doublemiss;
    } else {
      tmp[jobs] = static_cast<double>(vec[jobs]);
    }
  }
  this->putdb(col, tmp);
}

// -----------------------------------------------------------------------------

void ObsTable::putdb(const std::string & col, const std::vector<float> & vec) const {
  std::vector<double> tmp(vec.size());
  float floatmiss = util::missingValue(floatmiss);
  double doublemiss = util::missingValue(doublemiss);
  for (size_t jobs = 0; jobs < vec.size(); ++jobs) {
    if (vec[jobs] == floatmiss) {
      tmp[jobs] = doublemiss;
    } else {
      tmp[jobs] = static_cast<double>(vec[jobs]);
    }
  }
  this->putdb(col, tmp);
}

// -----------------------------------------------------------------------------

void ObsTable::putdb(const std::string & col, const std::vector<double> & vec) const {
  ASSERT(vec.size() == nobs());
  if (data_.find(col) != data_.end()) {
    oops::Log::info() << "ObsTable::putdb over-writing " << col << std::endl;
    data_[col] = vec;
  } else {
    data_.insert(std::pair<std::string, std::vector<double> >(col, vec));
  }
}

// -----------------------------------------------------------------------------

void ObsTable::getdb(const std::string & col, std::vector<int> & vec) const {
  std::vector<double> tmp;
  this->getdb(col, tmp);
  int intmiss = util::missingValue(intmiss);
  double doublemiss = util::missingValue(doublemiss);
  vec.resize(nobs());
  for (size_t jobs = 0; jobs < nobs(); ++jobs) {
    if (tmp[jobs] == doublemiss) {
      vec[jobs] = intmiss;
    } else {
      vec[jobs] = lround(tmp[jobs]);
    }
  }
}

// -----------------------------------------------------------------------------

void ObsTable::getdb(const std::string & col, std::vector<float> & vec) const {
  std::vector<double> tmp;
  this->getdb(col, tmp);
  float floatmiss = util::missingValue(floatmiss);
  double doublemiss = util::missingValue(doublemiss);
  vec.resize(nobs());
  for (size_t jobs = 0; jobs < nobs(); ++jobs) {
    if (tmp[jobs] == doublemiss) {
      vec[jobs] = floatmiss;
    } else {
      vec[jobs] = static_cast<float>(tmp[jobs]);
    }
  }
}

// -----------------------------------------------------------------------------

void ObsTable::getdb(const std::string & col, std::vector<double> & vec) const {
  std::map<std::string, std::vector<double> >::const_iterator ic = data_.find(col);
  if (ic == data_.end()) {
    oops::Log::error() << "ObsTable::getdb " << col << " not found." << std::endl;
    ABORT("ObsTable::getdb column not found");
  }
  vec.resize(nobs());
  for (unsigned int jobs = 0; jobs < nobs(); ++jobs) {
    vec[jobs] = ic->second[jobs];
  }
}

// -----------------------------------------------------------------------------

void ObsTable::generateDistribution(const eckit::Configuration & config) {
  oops::Log::trace() << "ObsTable::generateDistribution starting" << std::endl;

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

// Generate obs error
  const double err = config.getDouble("obs_error");
  std::vector<double> obserr(nobs);
  for (unsigned int jj = 0; jj < nobs; ++jj) {
    obserr[jj] = err;
  }
  this->putdb("ObsError", obserr);

  oops::Log::trace() << "ObsTable::generateDistribution done, nobs= " << nobs << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTable::random(std::vector<double> & data) const {
  util::NormalDistribution<double> x(data.size(), 0.0, 1.0, getSeed());
  for (size_t jj = 0; jj < data.size(); ++jj) data[jj] = x[jj];
}

// -----------------------------------------------------------------------------
//  ObsTable Private Methods
// -----------------------------------------------------------------------------

void ObsTable::otOpen(const std::string & filename) {
  oops::Log::trace() << "ObsTable::ot_read starting" << std::endl;
  std::ifstream fin(filename.c_str());
  if (!fin.is_open()) ABORT("ObsTable::otOpen: Error opening file: " + filename);

  int ncol, nobs;
  fin >> ncol;

  std::vector<std::string> colnames;
  for (int jc = 0; jc < ncol; ++jc) {
    std::string col;
    fin >> col;
    colnames.push_back(col);
  }

  fin >> nobs;
  locations_.clear();
  std::vector<double> newcol;
  for (int jc = 0; jc < ncol; ++jc) {
    ASSERT(data_.find(colnames[jc]) == data_.end());
    data_.insert(std::pair<std::string, std::vector<double> >(colnames[jc], newcol));
  }

  times_.clear();
  for (int jobs = 0; jobs < nobs; ++jobs) {
    int jjj;
    fin >> jjj;
    ASSERT(jjj == jobs);
    std::string sss;
    fin >> sss;
    util::DateTime ttt(sss);
    bool inside = ttt > winbgn_ && ttt <= winend_;

    if (inside) times_.push_back(ttt);
    double loc;
    fin >> loc;
    if (inside) locations_.push_back(loc);
    for (std::map<std::string, std::vector<double> >::iterator jo = data_.begin();
         jo != data_.end(); ++jo) {
      double val;
      fin >> val;
      if (inside) jo->second.push_back(val);
    }
  }

  fin.close();
  oops::Log::trace() << "ObsTable::ot_read done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTable::otWrite(const std::string & filename) const {
  oops::Log::trace() << "ObsTable::otWrite writing " << filename << std::endl;

  const size_t ioproc = 0;

  int nobs = times_.size();
  if (comm_.size() > 1) comm_.allReduceInPlace(nobs, eckit::mpi::Operation::SUM);

  std::vector<util::DateTime> timebuff(nobs);
  oops::mpi::gather(comm_, times_, timebuff, ioproc);

  std::vector<double> locbuff(nobs);
  oops::mpi::gather(comm_, locations_, locbuff, ioproc);

  std::vector<double> datasend(times_.size() * data_.size());
  size_t iobs = 0;
  for (size_t jobs = 0; jobs < times_.size(); ++jobs) {
    for (auto jo = data_.begin(); jo != data_.end(); ++jo) {
      datasend[iobs] = jo->second[jobs];
      ++iobs;
    }
  }
  std::vector<double> databuff(data_.size() * nobs);
  oops::mpi::gather(comm_, datasend, databuff, ioproc);

  if (comm_.rank() == ioproc) {
    std::ofstream fout(filename.c_str());
    if (!fout.is_open()) ABORT("ObsTable::otWrite: Error opening file: " + filename);

    int ncol = data_.size();
    fout << ncol << std::endl;

    for (auto jo = data_.begin(); jo != data_.end(); ++jo)
      fout << jo->first << std::endl;

    fout << nobs << std::endl;

    size_t iii = 0;
    for (int jobs = 0; jobs < nobs; ++jobs) {
      fout << jobs;
      fout << "  " << timebuff[jobs];
      fout << "  " << locbuff[jobs];
      for (int jcol = 0; jcol < ncol; ++jcol) {
        fout << "  " << databuff[iii];
        ++iii;
      }
      fout << std::endl;
    }

    fout.close();
  }

  oops::Log::trace() << "ObsTable::otWrite done" << std::endl;
}

// -----------------------------------------------------------------------------
ObsIterator ObsTable::begin() const {
  return ObsIterator(locations_, 0);
}
// -----------------------------------------------------------------------------
ObsIterator ObsTable::end() const {
  return ObsIterator(locations_, locations_.size());
}
// -----------------------------------------------------------------------------

void ObsTable::print(std::ostream & os) const {
  os << "ObsTable: assimilation window = " << winbgn_ << " to " << winend_ << std::endl;
  os << "ObsTable: file in = " << nameIn_ << ", file out = " << nameOut_;
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
