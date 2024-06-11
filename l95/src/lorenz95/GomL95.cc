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

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iomanip>
#include <limits>

#include "eckit/config/Configuration.h"

#include "lorenz95/L95Traits.h"
#include "lorenz95/LocsL95.h"
#include "lorenz95/ObsTable.h"
#include "oops/base/Locations.h"
#include "oops/base/Variable.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/Logger.h"
#include "oops/util/Random.h"

namespace oops {
class Variables;
}

namespace lorenz95 {

// -----------------------------------------------------------------------------
GomL95::GomL95(const oops::Locations<L95ObsTraits> & locs,
               const oops::Variables &, const std::vector<size_t> &)
  : size_(locs.samplingMethod(oops::Variable("x")).sampledLocations().size()), locval_(size_, 0.0)
{
  oops::Log::trace() << "GomL95::GomL95 done" << std::endl;
}
// -----------------------------------------------------------------------------
/*! Constructor with Configuration */
GomL95::GomL95(const eckit::Configuration & conf, const ObsTable &, const oops::Variables &)
  : size_(0), locval_()
{
  this->read(conf);
}
// -----------------------------------------------------------------------------
GomL95 & GomL95::operator*=(const double & zz) {
  for (size_t jj = 0; jj < size_; ++jj) locval_[jj] *= zz;
  return *this;
}
// -----------------------------------------------------------------------------
GomL95 & GomL95::operator+=(const GomL95 & rhs)
{
  for (size_t jj = 0; jj < size_; ++jj) locval_[jj] += rhs.locval_[jj];
  return *this;
}
// -----------------------------------------------------------------------------
GomL95 & GomL95::operator-=(const GomL95 & rhs)
{
  for (size_t jj = 0; jj < size_; ++jj) locval_[jj] -= rhs.locval_[jj];
  return *this;
}
// -----------------------------------------------------------------------------
GomL95 & GomL95::operator*=(const GomL95 & rhs)
{
  for (size_t jj = 0; jj < size_; ++jj) locval_[jj] *= rhs.locval_[jj];
  return *this;
}
// -----------------------------------------------------------------------------
void GomL95::zero() {
  for (size_t jj = 0; jj < size_; ++jj) locval_[jj] = 0.0;
}
// -----------------------------------------------------------------------------
double GomL95::rms() const {
  double xnorm(0.0);
  for (size_t jj = 0; jj < size_; ++jj) xnorm += locval_[jj] * locval_[jj];
  return sqrt(xnorm/static_cast<double>(size_));
}
// -----------------------------------------------------------------------------
double GomL95::normalizedrms(const GomL95 & rhs) const {
  GomL95 temp_gv(*this);
  for (size_t jj = 0; jj < size_; ++jj) temp_gv.locval_[jj] /= rhs.locval_[jj];
  return temp_gv.rms();
}
// -----------------------------------------------------------------------------
void GomL95::random() {
  util::NormalDistribution<double> x(size_, 0.0, 1.0, 5);
  for (size_t jj = 0; jj < size_; ++jj) locval_[jj] = x[jj];
}
// -----------------------------------------------------------------------------
double GomL95::dot_product_with(const GomL95 & gom) const {
  double zz = 0.0;
  for (size_t jj = 0; jj < size_; ++jj) zz += locval_[jj] * gom.locval_[jj];
  return zz;
}
// -----------------------------------------------------------------------------
void GomL95::fill(const oops::Variable &, const ConstVectorRef<size_t> &indx,
                  const ConstMatrixRef<double> &vals, const bool) {
  ASSERT(indx.size() == vals.rows());
  ASSERT(vals.cols() == 1);
  for (Eigen::Index jj = 0; jj < indx.size(); ++jj) locval_[indx[jj]] = vals(jj, 0);
}
// -----------------------------------------------------------------------------
void GomL95::fillAD(const oops::Variable &, const ConstVectorRef<size_t> &indx,
                    MatrixRef<double> vals, const bool) const {
  ASSERT(indx.size() == vals.rows());
  ASSERT(vals.cols() == 1);
  for (Eigen::Index jj = 0; jj < indx.size(); ++jj) vals(jj, 0) += locval_[indx[jj]];
}
// -----------------------------------------------------------------------------
void GomL95::read(const eckit::Configuration & conf) {
  const std::string filename(conf.getString("filename"));
  oops::Log::trace() << "GomL95::read opening " << filename << std::endl;
  std::ifstream fin(filename.c_str());
  if (!fin.is_open()) ABORT("GomL95::read: Error opening file: " + filename);

  size_t size;
  fin >> size;

  if (size_ != size) {
    size_ = size;
    locval_.resize(size_);
  }

  for (size_t jj = 0; jj < size_; ++jj) fin >> locval_[jj];

  fin.close();
  oops::Log::trace() << "GomL95::read: file closed." << std::endl;
}
// -----------------------------------------------------------------------------
void GomL95::write(const eckit::Configuration & conf) const {
  const std::string filename(conf.getString("filename"));
  oops::Log::trace() << "GomL95::write opening " << filename << std::endl;
  std::ofstream fout(filename.c_str());
  if (!fout.is_open()) ABORT("GomL95::write: Error opening file: " + filename);

  fout << size_ << std::endl;
  fout << std::endl;
  fout.precision(std::numeric_limits<double>::digits10);
  for (size_t jj = 0; jj < size_; ++jj) fout << locval_[jj] << " ";
  fout << std::endl;

  fout.close();
  oops::Log::trace() << "GomL95::write file closed." << std::endl;
}
// -----------------------------------------------------------------------------
void GomL95::print(std::ostream & os) const {
  if (size_ > 0) {
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
  } else {
    os << " No observations";
  }
}
// -----------------------------------------------------------------------------

}  // namespace lorenz95
