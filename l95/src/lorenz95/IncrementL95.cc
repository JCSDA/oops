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

#include "lorenz95/IncrementL95.h"

#include <fstream>
#include <string>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/LocalIncrement.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/stringFunctions.h"

#include "lorenz95/FieldL95.h"
#include "lorenz95/GomL95.h"
#include "lorenz95/Iterator.h"
#include "lorenz95/ModelBiasCorrection.h"
#include "lorenz95/Resolution.h"
#include "lorenz95/StateL95.h"

namespace oops {
  class Variables;
}
namespace sf = util::stringfunctions;

namespace lorenz95 {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
IncrementL95::IncrementL95(const Resolution & resol, const oops::Variables & vars,
                           const util::DateTime & vt)
  : fld_(resol), time_(vt), vars_(vars)
{
  fld_.zero();
  oops::Log::trace() << "IncrementL95::IncrementL95 created." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementL95::IncrementL95(const Resolution & resol, const IncrementL95 & dx)
  : fld_(resol), time_(dx.time_), vars_(dx.variables())
{
  fld_ = dx.fld_;
  oops::Log::trace() << "IncrementL95::IncrementL95 created by interpolation." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementL95::IncrementL95(const IncrementL95 & dx, const bool copy)
  : fld_(dx.fld_, copy), time_(dx.time_), vars_(dx.variables())
{
  oops::Log::trace() << "IncrementL95::IncrementL95 copy-created." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementL95::~IncrementL95() {
  oops::Log::trace() << "IncrementL95::~IncrementL95 destructed" << std::endl;
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
void IncrementL95::ones() {
  fld_.ones();
}
// -----------------------------------------------------------------------------
void IncrementL95::dirac(const eckit::Configuration & conf) {
  fld_.dirac(conf);
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
void IncrementL95::random(const size_t & seed) {
  fld_.random(seed);
}
// -----------------------------------------------------------------------------
void IncrementL95::accumul(const double & zz, const StateL95 & xx) {
  fld_.axpy(zz, xx.getField());
}
// -----------------------------------------------------------------------------
/// Utilities
// -----------------------------------------------------------------------------
void IncrementL95::read(const eckit::Configuration & config) {
  std::string filename(config.getString("filename"));
  sf::swapNameMember(config, filename);
  oops::Log::trace() << "IncrementL95::read opening " << filename << std::endl;
  std::ifstream fin(filename.c_str());
  if (!fin.is_open()) ABORT("IncrementL95::read: Error opening file: " + filename);

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
  oops::Log::trace() << "IncrementL95::read: file closed." << std::endl;
}
// -----------------------------------------------------------------------------
void IncrementL95::write(const eckit::Configuration & config) const {
  const std::string dir = config.getString("datadir");
  const std::string exp = config.getString("exp");
  const std::string type = config.getString("type");
  std::string filename = dir+"/"+exp+"."+type;

  if (type == "krylov") {
    if (!config.has("iteration"))
      throw eckit::BadValue("'iteration' was not set in the parameters passed to write() "
                            "even though 'type' was set to '" + type + "'", Here());
    const int iter = config.getInt("iteration");
    filename += "."+std::to_string(iter);
  }

  const util::DateTime antime(config.getString("date"));
  filename += "."+antime.toString();
  const util::Duration step = time_ - antime;
  filename += "."+step.toString();
  sf::swapNameMember(config, filename);
  filename += ".l95";

  oops::Log::trace() << "IncrementL95::write opening " << filename << std::endl;
  std::ofstream fout(filename.c_str());
  if (!fout.is_open()) ABORT("IncrementL95::write: Error opening file: " + filename);

  fout << fld_.resol() << std::endl;
  fout << time_ << std::endl;
  fld_.write(fout);
  fout << std::endl;

  fout.close();
  oops::Log::trace() << "IncrementL95::write file closed." << std::endl;
}
// -----------------------------------------------------------------------------
void IncrementL95::print(std::ostream & os) const {
  os << std::endl << " Valid time: " << time_;
  os << std::endl << fld_;
}
// -----------------------------------------------------------------------------
oops::LocalIncrement IncrementL95::getLocal(const Iterator & i) const {
  std::vector<std::string> vars;
  vars.push_back("x");
  std::vector<double> vals;
  vals.push_back(fld_[i.index()]);
  std::vector<int> varlens;
  varlens.push_back(1);
  return oops::LocalIncrement(oops::Variables(vars), vals, varlens);
}
// -----------------------------------------------------------------------------
void IncrementL95::setLocal(const oops::LocalIncrement & gp, const Iterator & i) {
  std::vector<double> vals;
  vals = gp.getVals();
  fld_[i.index()] = vals[0];
}
// -----------------------------------------------------------------------------
/// Serialize - deserialize
// -----------------------------------------------------------------------------
size_t IncrementL95::serialSize() const {
  size_t nn = 0;
  nn += fld_.serialSize();
  nn += time_.serialSize();
  return nn;
}
// -----------------------------------------------------------------------------
void IncrementL95::serialize(std::vector<double> & vect) const {
  fld_.serialize(vect);
  time_.serialize(vect);
}
// -----------------------------------------------------------------------------
void IncrementL95::deserialize(const std::vector<double> & vect, size_t & index) {
  const size_t ii = index + this->serialSize();
  fld_.deserialize(vect, index);
  time_.deserialize(vect, index);
  ASSERT(index == ii);
}
// -----------------------------------------------------------------------------

}  // namespace lorenz95
