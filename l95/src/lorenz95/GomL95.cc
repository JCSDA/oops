/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/GomL95.h"

#include <cstdlib>
#include <limits>
#include <fstream>
#include <random>
#include <vector>

#include "eckit/config/Configuration.h"
#include "lorenz95/LocsL95.h"
#include "lorenz95/ObsTable.h"
#include "util/abor1_cpp.h"
#include "util/Logger.h"

namespace oops {
  class Variables;
}

namespace lorenz95 {

// -----------------------------------------------------------------------------
GomL95::GomL95(const LocsL95 & locs, const oops::Variables &)
  : size_(0), iobs_(), locval_(), current_(0)
{
  oops::Log::trace() << "GomL95::GomL95 starting " << std::endl;
  size_ = locs.size();
  iobs_.resize(size_);
  locval_.resize(size_);
  for (size_t jj = 0; jj < size_; ++jj) iobs_[jj] = locs.globalIndex(jj);
  for (size_t jj = 0; jj < size_; ++jj) locval_[jj] = locs[jj];
}
// -----------------------------------------------------------------------------
GomL95::GomL95(const eckit::Configuration & conf, const oops::Variables &)
  : size_(0), iobs_(), locval_(), current_(0)
{
  this->read(conf);
}
// -----------------------------------------------------------------------------
GomL95::~GomL95() {}
// -----------------------------------------------------------------------------
GomL95 & GomL95::operator*=(const double & zz) {
  for (size_t jj = 0; jj < size_; ++jj) locval_[jj] *= zz;
  return *this;
}
// -----------------------------------------------------------------------------
void GomL95::zero() {
  for (size_t jj = 0; jj < size_; ++jj) locval_[jj] = 0.0;
}
// -----------------------------------------------------------------------------
void GomL95::random() {
  static std::mt19937 generator(5);
  static std::normal_distribution<double> distribution(0.0, 1.0);
  for (size_t jj = 0; jj < size_; ++jj) locval_[jj] = distribution(generator);
}
// -----------------------------------------------------------------------------
double GomL95::dot_product_with(const GomL95 & gom) const {
  double zz = 0.0;
  for (size_t jj = 0; jj < size_; ++jj) zz += locval_[jj] * gom.locval_[jj];
  return zz;
}
// -----------------------------------------------------------------------------
void GomL95::read(const eckit::Configuration & conf) {
  const std::string filename(conf.getString("filename"));
  oops::Log::trace() << "GomL95::read opening " << filename << std::endl;
  std::ifstream fin(filename.c_str());
  if (!fin.is_open()) ABORT("GomL95::read: Error opening file");

  size_t size;
  fin >> size;

  if (size_ != size) {
    size_ = size;
    iobs_.resize(size_);
    locval_.resize(size_);
  }

  for (size_t jj = 0; jj < size_; ++jj) fin >> iobs_[jj];
  for (size_t jj = 0; jj < size_; ++jj) fin >> locval_[jj];

  fin.close();
  oops::Log::trace() << "GomL95::read: file closed." << std::endl;
}
// -----------------------------------------------------------------------------
void GomL95::write(const eckit::Configuration & conf) const {
  const std::string filename(conf.getString("filename"));
  oops::Log::trace() << "GomL95::write opening " << filename << std::endl;
  std::ofstream fout(filename.c_str());
  if (!fout.is_open()) ABORT("GomL95::write: Error opening file");

  fout << size_ << std::endl;
  for (size_t jj = 0; jj < size_; ++jj) fout << iobs_[jj] << " ";
  fout << std::endl;
  fout.precision(std::numeric_limits<double>::digits10);
  for (size_t jj = 0; jj < size_; ++jj) fout << locval_[jj] << " ";
  fout << std::endl;

  fout.close();
  oops::Log::trace() << "GomL95::write file closed." << std::endl;
}
// -----------------------------------------------------------------------------
void GomL95::print(std::ostream & os) const {
  double zmin = locval_[0];
  double zmax = locval_[0];
  double zavg = 0.0;
  for (size_t jj = 0; jj < size_; ++jj) {
    if (locval_[jj] < zmin) zmin = locval_[jj];
    if (locval_[jj] > zmax) zmax = locval_[jj];
    zavg += locval_[jj];
  }
  zavg /= size_;
  os << size_ << "values,  Min=" << zmin << ", Max=" << zmax << ", Average=" << zavg;
}
// -----------------------------------------------------------------------------

}  // namespace lorenz95
