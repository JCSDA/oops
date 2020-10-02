/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#include "lorenz95/ObsTableView.h"

#include <algorithm>
#include <cmath>
#include <iostream>
#include <numeric>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/types/Types.h"

#include "oops/util/Logger.h"
#include "oops/util/missingValues.h"

// -----------------------------------------------------------------------------
namespace lorenz95 {

ObsTableView::ObsTableView(const eckit::Configuration & config, const eckit::mpi::Comm & comm,
                           const util::DateTime & bgn, const util::DateTime & end)
  : obstable_(new ObsTable(config, comm, bgn, end)),
    localobs_(obstable_->nobs()), obsdist_(obstable_->nobs(), 0.0)
{
  std::iota(localobs_.begin(), localobs_.end(), 0);
  oops::Log::trace() << "ObsTableView::ObsTableView created nobs = " << nobs() << std::endl;
}

// -----------------------------------------------------------------------------

ObsTableView::ObsTableView(const ObsTableView & obstable,
                           const eckit::geometry::Point2 & center,
                           const eckit::Configuration & conf)
  :  obstable_(obstable.obstable_), localobs_(), obsdist_()
{
  std::vector<double> locations = obstable.locations();
  const double dist = conf.getDouble("lengthscale");
  for (unsigned int jj = 0; jj < obstable.nobs(); ++jj) {
    double curdist = std::abs(center[0] - locations[jj]);
    curdist = std::min(curdist, 1.-curdist);
    if ( curdist < dist ) {
      localobs_.push_back(jj);
      obsdist_.push_back(curdist);
    }
  }
  oops::Log::trace() << "ObsTableView::ObsTableView created" << std::endl;
}

// -----------------------------------------------------------------------------

ObsTableView::~ObsTableView() {
  oops::Log::trace() << "ObsTableView::ObsTableView destructed" << std::endl;
}

// -----------------------------------------------------------------------------

bool ObsTableView::has(const std::string & col) const {
  oops::Log::trace() << "ObsTableView::has" << std::endl;
  return obstable_->has(col);
}

// -----------------------------------------------------------------------------

void ObsTableView::putdb(const std::string & col, const std::vector<int> & vec) const {
  int missing;
  std::vector<int> fullvec(obstable_->nobs(), util::missingValue(missing));
  if (obstable_->has(col)) {
    obstable_->getdb(col, fullvec);
  }
  for (unsigned int i = 0; i < nobs(); i++) {
    fullvec[localobs_[i]] = vec[i];
  }
  obstable_->putdb(col, fullvec);
  oops::Log::trace() << "ObsTableView::putdb done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTableView::putdb(const std::string & col, const std::vector<float> & vec) const {
  float missing;
  std::vector<float> fullvec(obstable_->nobs(), util::missingValue(missing));
  if (obstable_->has(col)) {
    obstable_->getdb(col, fullvec);
  }
  for (unsigned int i = 0; i < nobs(); i++) {
    fullvec[localobs_[i]] = vec[i];
  }
  obstable_->putdb(col, fullvec);
  oops::Log::trace() << "ObsTableView::putdb done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTableView::putdb(const std::string & col, const std::vector<double> & vec) const {
  double missing;
  std::vector<double> fullvec(obstable_->nobs(), util::missingValue(missing));
  if (obstable_->has(col)) {
    obstable_->getdb(col, fullvec);
  }
  for (unsigned int i = 0; i < nobs(); i++) {
    fullvec[localobs_[i]] = vec[i];
  }
  obstable_->putdb(col, fullvec);
  oops::Log::trace() << "ObsTableView::putdb done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTableView::getdb(const std::string & col, std::vector<int> & vec) const {
  std::vector<int> fullvec;
  obstable_->getdb(col, fullvec);
  vec.resize(nobs());
  for (unsigned int i = 0; i < nobs(); i++) {
    vec[i] = fullvec[localobs_[i]];
  }
  oops::Log::trace() << "ObsTableView::getdb done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTableView::getdb(const std::string & col, std::vector<float> & vec) const {
  std::vector<float> fullvec;
  obstable_->getdb(col, fullvec);
  vec.resize(nobs());
  for (unsigned int i = 0; i < nobs(); i++) {
    vec[i] = fullvec[localobs_[i]];
  }
  oops::Log::trace() << "ObsTableView::getdb done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTableView::getdb(const std::string & col, std::vector<double> & vec) const {
  std::vector<double> fullvec;
  obstable_->getdb(col, fullvec);
  vec.resize(nobs());
  for (unsigned int i = 0; i < nobs(); i++) {
    vec[i] = fullvec[localobs_[i]];
  }
  oops::Log::trace() << "ObsTableView::getdb done" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTableView::random(std::vector<double> & v) const {
  obstable_->random(v);
  oops::Log::trace() << "ObsTableView::random done" << std::endl;
}

// -----------------------------------------------------------------------------

unsigned int ObsTableView::nobs() const {
  return localobs_.size();
  oops::Log::trace() << "ObsTableView::nobs done" << std::endl;
}

// -----------------------------------------------------------------------------

std::vector<double> ObsTableView::locations() const {
  std::vector<double> full = obstable_->locations();
  std::vector<double> local(nobs());
  for (unsigned int i = 0; i < nobs(); i++) {
    local[i] = full[localobs_[i]];
  }
  oops::Log::trace() << "ObsTableView::locations done" << std::endl;
  return local;
}

// -----------------------------------------------------------------------------

void ObsTableView::generateDistribution(const eckit::Configuration & conf) {
  obstable_->generateDistribution(conf);
  int nobs = obstable_->nobs();
  for (int i = 0; i < nobs; i++)
    localobs_.push_back(i);
  oops::Log::trace() << "ObsTableView::generateDistribution done" << std::endl;
}

// -----------------------------------------------------------------------------

std::unique_ptr<LocsL95> ObsTableView::locations(const util::DateTime & t1,
                         const util::DateTime & t2) const {
  // get times and locations from the obsspace
  std::vector<util::DateTime> all_times = obstable_->times();
  std::vector<double> all_locs = obstable_->locations();
  // find local times that are within t1 and t2
  std::vector<int> mask;
  for (unsigned int i = 0; i < nobs(); i++) {
    if (all_times[localobs_[i]] > t1 && all_times[localobs_[i]] <= t2)
      mask.push_back(i);
  }
  // set up locations
  const unsigned int nobs_t = mask.size();
  std::vector<double> locs(nobs_t);
  std::vector<util::DateTime> times(nobs_t);
  for (unsigned int i = 0; i < nobs_t; i++) {
    locs[i] = all_locs[localobs_[mask[i]]];
    times[i] = all_times[localobs_[mask[i]]];
  }
  oops::Log::trace() << "ObsTableView::locations done" << std::endl;
  return std::unique_ptr<LocsL95>(new LocsL95(locs, times));
}

// -----------------------------------------------------------------------------

void ObsTableView::printJo(const ObsVec1D & x1, const ObsVec1D & x2) {
  oops::Log::info() << "ObsTableView::printJo not implemented" << std::endl;
}

// -----------------------------------------------------------------------------

void ObsTableView::print(std::ostream & os) const {
  os << "Local observation indices: " << localobs_ << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace lorenz95
