/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "lorenz95/IncrementL95.h"

#include <fstream>
#include <string>

#include "util/Logger.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/abor1_cpp.h"
#include "util/dot_product.h"

#include "lorenz95/FieldL95.h"
#include "lorenz95/GomL95.h"
#include "lorenz95/LocsL95.h"
#include "lorenz95/ModelBiasCorrection.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/StateL95.h"

using oops::Log;

namespace lorenz95 {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
IncrementL95::IncrementL95(const Resolution & resol, const NoVariables &,
                           const util::DateTime & vt)
  : fld_(resol), time_(vt)
{
  fld_.zero();
  Log::trace() << "IncrementL95::IncrementL95 created." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementL95::IncrementL95(const Resolution & resol, const IncrementL95 & dx)
  : fld_(resol), time_(dx.time_)
{
  fld_ = dx.fld_;
  Log::trace() << "IncrementL95::IncrementL95 created by interpolation." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementL95::IncrementL95(const IncrementL95 & dx, const bool copy)
  : fld_(dx.fld_), time_(dx.time_)
{
  Log::trace() << "IncrementL95::IncrementL95 copy-created." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementL95::~IncrementL95() {
  Log::trace() << "IncrementL95::~IncrementL95 destructed" << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
void IncrementL95::diff(const StateL95 & x1, const StateL95 & x2) {
  ASSERT(time_ == x1.validTime());
  ASSERT(time_ == x2.validTime());
  fld_.diff(x1.getField(), x2.getField());
}
// -----------------------------------------------------------------------------
IncrementL95 & IncrementL95::operator=(const IncrementL95 & rhs) {
  fld_ = rhs.fld_;
  time_ = rhs.time_;
  return *this;
}
// -----------------------------------------------------------------------------
IncrementL95 & IncrementL95::operator+=(const IncrementL95 & rhs) {
  ASSERT(time_ == rhs.time_);
  fld_ += rhs.fld_;
  return *this;
}
// -----------------------------------------------------------------------------
IncrementL95 & IncrementL95::operator-=(const IncrementL95 & rhs) {
  ASSERT(time_ == rhs.time_);
  fld_ -= rhs.fld_;
  return *this;
}
// -----------------------------------------------------------------------------
IncrementL95 & IncrementL95::operator*=(const double & zz) {
  fld_ *= zz;
  return *this;
}
// -----------------------------------------------------------------------------
void IncrementL95::zero() {
  fld_.zero();
}
// -----------------------------------------------------------------------------
void IncrementL95::zero(const util::DateTime & vt) {
  fld_.zero();
  time_ = vt;
}
// -----------------------------------------------------------------------------
void IncrementL95::axpy(const double & zz, const IncrementL95 & rhs,
                        const bool check) {
  ASSERT(!check || time_ == rhs.time_);
  fld_.axpy(zz, rhs.fld_);
}
// -----------------------------------------------------------------------------
double IncrementL95::dot_product_with(const IncrementL95 & other) const {
  double zz = dot_product(fld_, other.fld_);
  return zz;
}
// -----------------------------------------------------------------------------
void IncrementL95::schur_product_with(const IncrementL95 & rhs) {
  fld_.schur(rhs.fld_);
}
// -----------------------------------------------------------------------------
void IncrementL95::random() {
  fld_.random();
}
// -----------------------------------------------------------------------------
void IncrementL95::accumul(const double & zz, const StateL95 & xx) {
  fld_.axpy(zz, xx.getField());
}
// -----------------------------------------------------------------------------
/// Utilities
// -----------------------------------------------------------------------------
void IncrementL95::read(const eckit::Configuration & config) {
  const std::string filename(config.getString("filename"));
  Log::trace() << "IncrementL95::read opening " << filename << std::endl;
  std::ifstream fin(filename.c_str());
  if (!fin.is_open()) ABORT("IncrementL95::read: Error opening file");

  int resol;
  fin >> resol;
  ASSERT(fld_.resol() == resol);

  std::string stime;
  fin >> stime;
  const util::DateTime tt(stime);
  const util::DateTime tc(config.getString("date"));
  if (tc != tt) {
    ABORT("IncrementL95::read: date and data file inconsistent.");
  }
  time_ = tt;

  fld_.read(fin);

  fin.close();
  Log::trace() << "IncrementL95::read: file closed." << std::endl;
}
// -----------------------------------------------------------------------------
void IncrementL95::write(const eckit::Configuration & config) const {
  std::string dir = config.getString("datadir");
  std::string exp = config.getString("exp");
  std::string type = config.getString("type");
  std::string filename = dir+"/"+exp+"."+type;

  const util::DateTime antime(config.getString("date"));
  filename += "."+antime.toString();
  const util::Duration step = time_ - antime;
  filename += "."+step.toString();

  Log::trace() << "IncrementL95::write opening " << filename << std::endl;
  std::ofstream fout(filename.c_str());
  if (!fout.is_open()) ABORT("IncrementL95::write: Error opening file");

  fout << fld_.resol() << std::endl;
  fout << time_ << std::endl;
  fld_.write(fout);
  fout << std::endl;

  fout.close();
  Log::trace() << "IncrementL95::write file closed." << std::endl;
}
// -----------------------------------------------------------------------------
void IncrementL95::print(std::ostream & os) const {
  os << std::endl << " Valid time: " << time_;
  os << std::endl << fld_;
}
// -----------------------------------------------------------------------------
/// Interpolate to observation location
// -----------------------------------------------------------------------------
void IncrementL95::interpolateTL(const LocsL95 & locs, GomL95 & vals) const {
  fld_.interp(locs, vals);
}
// -----------------------------------------------------------------------------
void IncrementL95::interpolateAD(const LocsL95 & locs, const GomL95 & vals) {
  fld_.interpAD(locs, vals);
}
// -----------------------------------------------------------------------------

}  // namespace lorenz95
