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

#include "atlas/field.h"

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
#include "lorenz95/LocsL95.h"
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
IncrementL95::IncrementL95(const Resolution & resol, const oops::Variables &,
                           const util::DateTime & vt)
  : fld_(resol), time_(vt)
{
  fld_.zero();
  oops::Log::trace() << "IncrementL95::IncrementL95 created." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementL95::IncrementL95(const Resolution & resol, const IncrementL95 & dx)
  : fld_(resol), time_(dx.time_)
{
  fld_ = dx.fld_;
  oops::Log::trace() << "IncrementL95::IncrementL95 created by interpolation." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementL95::IncrementL95(const IncrementL95 & dx, const bool copy)
  : fld_(dx.fld_, copy), time_(dx.time_)
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
void IncrementL95::dirac(const DiracParameters_ & params) {
  fld_.dirac(params);
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
void IncrementL95::read(const ReadParameters_ & params) {
  std::string filename(params.filename);
  sf::swapNameMember(params.member, filename);
  oops::Log::trace() << "IncrementL95::read opening " << filename << std::endl;
  std::ifstream fin(filename.c_str());
  if (!fin.is_open()) ABORT("IncrementL95::read: Error opening file: " + filename);

  int resol;
  fin >> resol;
  ASSERT(fld_.resol() == resol);

  std::string stime;
  fin >> stime;
  const util::DateTime tt(stime);
  const util::DateTime tc(params.date);
  if (tc != tt) {
    ABORT("IncrementL95::read: date and data file inconsistent.");
  }
  time_ = tt;

  fld_.read(fin);

  fin.close();
  oops::Log::trace() << "IncrementL95::read: file closed." << std::endl;
}
// -----------------------------------------------------------------------------
void IncrementL95::write(const WriteParameters_ & params) const {
  const std::string &dir = params.datadir;
  const std::string &exp = params.exp;
  const std::string &type = params.type;
  std::string filename = dir+"/"+exp+"."+type;

  if (type == "krylov") {
    if (params.iteration.value() == boost::none)
      throw eckit::BadValue("'iteration' was not set in the parameters passed to write() "
                            "even though 'type' was set to '" + type + "'", Here());
    const int &iter = *params.iteration.value();
    filename += "."+std::to_string(iter);
  }

  if (params.date.value() == boost::none)
    throw eckit::BadValue("'date' was not set in the parameters passed to write()", Here());
  const util::DateTime &antime = *params.date.value();
  filename += "."+antime.toString();
  const util::Duration step = time_ - antime;
  filename += "."+step.toString();
  sf::swapNameMember(params.member, filename);
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
/// Convert to/from ATLAS fieldset
// -----------------------------------------------------------------------------
void IncrementL95::setAtlas(atlas::FieldSet *) const {
  ABORT("FieldL95 setAtlas not implemented");
}
// -----------------------------------------------------------------------------
void IncrementL95::toAtlas(atlas::FieldSet *) const {
  ABORT("FieldL95 toAtlas not implemented");
}
// -----------------------------------------------------------------------------
void IncrementL95::fromAtlas(atlas::FieldSet *) {
  ABORT("FieldL95 fromAtlas not implemented");
}
// -----------------------------------------------------------------------------
/// Serialize - deserialize
// -----------------------------------------------------------------------------
size_t IncrementL95::serialSize() const {
  size_t nn = 3;
  nn += fld_.serialSize();
  nn += time_.serialSize();
  return nn;
}
// -----------------------------------------------------------------------------
void IncrementL95::serialize(std::vector<double> & vect) const {
  vect.push_back(1000.0);
  fld_.serialize(vect);
  vect.push_back(2000.0);
  time_.serialize(vect);
  vect.push_back(3000.0);
}
// -----------------------------------------------------------------------------
void IncrementL95::deserialize(const std::vector<double> & vect, size_t & index) {
  size_t ii = index + this->serialSize();
  ASSERT(vect.at(index) == 1000.0);
  ++index;
  fld_.deserialize(vect, index);
  ASSERT(vect.at(index) == 2000.0);
  ++index;
  time_.deserialize(vect, index);
  ASSERT(vect.at(index) == 3000.0);
  ++index;
  ASSERT(index == ii);
}
// -----------------------------------------------------------------------------

}  // namespace lorenz95
