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

#include "model/IncrementQG.h"

#include <algorithm>
#include <string>
#include <utility>
#include <vector>

#include "atlas/field.h"

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "model/ErrorCovarianceQG.h"
#include "model/FieldsQG.h"
#include "model/GeometryQG.h"
#include "model/GomQG.h"
#include "model/ModelBiasIncrement.h"
#include "model/StateQG.h"

namespace qg {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
IncrementQG::IncrementQG(const GeometryQG & resol, const oops::Variables & vars,
                         const util::DateTime & vt)
  : fields_(new FieldsQG(resol, vars, lbc_, vt))
{
  fields_->zero();
  oops::Log::trace() << "IncrementQG constructed." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementQG::IncrementQG(const GeometryQG & resol, const IncrementQG & other, const bool ad)
  : fields_(new FieldsQG(*other.fields_, resol, ad))
{
  oops::Log::trace() << "IncrementQG constructed from other." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementQG::IncrementQG(const IncrementQG & other, const bool copy)
  : fields_(new FieldsQG(*other.fields_, copy))
{
  oops::Log::trace() << "IncrementQG copy-created." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementQG::IncrementQG(const IncrementQG & other)
  : fields_(new FieldsQG(*other.fields_))
{
  oops::Log::trace() << "IncrementQG copy-created." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementQG::~IncrementQG() {
  oops::Log::trace() << "IncrementQG destructed" << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
void IncrementQG::diff(const StateQG & x1, const StateQG & x2) {
  ASSERT(this->validTime() == x1.validTime());
  ASSERT(this->validTime() == x2.validTime());
  fields_->diff(x1.fields(), x2.fields());
}
// -----------------------------------------------------------------------------
IncrementQG & IncrementQG::operator=(const IncrementQG & rhs) {
  *fields_ = *rhs.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
IncrementQG & IncrementQG::operator+=(const IncrementQG & dx) {
  ASSERT(this->validTime() == dx.validTime());
  *fields_ += *dx.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
IncrementQG & IncrementQG::operator-=(const IncrementQG & dx) {
  ASSERT(this->validTime() == dx.validTime());
  *fields_ -= *dx.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
IncrementQG & IncrementQG::operator*=(const double & zz) {
  *fields_ *= zz;
  return *this;
}
// -----------------------------------------------------------------------------
void IncrementQG::zero() {
  fields_->zero();
}
// -----------------------------------------------------------------------------
void IncrementQG::zero(const util::DateTime & vt) {
  fields_->zero(vt);
}
// -----------------------------------------------------------------------------
void IncrementQG::ones() {
  fields_->ones();
}
// -----------------------------------------------------------------------------
void IncrementQG::axpy(const double & zz, const IncrementQG & dx,
                       const bool check) {
  ASSERT(!check || this->validTime() == dx.validTime());
  fields_->axpy(zz, *dx.fields_);
}
// -----------------------------------------------------------------------------
void IncrementQG::accumul(const double & zz, const StateQG & xx) {
  fields_->axpy(zz, xx.fields());
}
// -----------------------------------------------------------------------------
void IncrementQG::schur_product_with(const IncrementQG & dx) {
  fields_->schur_product_with(*dx.fields_);
}
// -----------------------------------------------------------------------------
double IncrementQG::dot_product_with(const IncrementQG & other) const {
  return dot_product(*fields_, *other.fields_);
}
// -----------------------------------------------------------------------------
void IncrementQG::random() {
  fields_->random();
}
// -----------------------------------------------------------------------------
void IncrementQG::dirac(const eckit::Configuration & config) {
  fields_->zero();
  if (config.has("date")) {
    util::DateTime dd(config.getString("date"));
    if (this->validTime() == dd) fields_->dirac(config);
  } else if (config.has("dates")) {
    std::vector<std::string> dates(config.getStringVector("dates"));
    for (const auto & date : dates) {
      util::DateTime dd(date);
      if (this->validTime() == dd) fields_->dirac(config);
    }
  } else {
    ABORT("missing date or dates in dirac configuration");
  }
}
// -----------------------------------------------------------------------------
/// ATLAS FieldSet
// -----------------------------------------------------------------------------
void IncrementQG::toFieldSet(atlas::FieldSet & fset) const {
  fields_->toFieldSet(fset);
}
// -----------------------------------------------------------------------------
void IncrementQG::fromFieldSet(const atlas::FieldSet & fset) {
  fields_->fromFieldSet(fset);
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void IncrementQG::read(const eckit::Configuration & files) {
  fields_->read(files);
}
// -----------------------------------------------------------------------------
void IncrementQG::write(const eckit::Configuration & files) const {
  fields_->write(files);
}
// -----------------------------------------------------------------------------
/// Serialization
// -----------------------------------------------------------------------------
size_t IncrementQG::serialSize() const {
  size_t nn = fields_->serialSize();
  return nn;
}
// -----------------------------------------------------------------------------
void IncrementQG::serialize(std::vector<double> & vect) const {
  fields_->serialize(vect);
}
// -----------------------------------------------------------------------------
void IncrementQG::deserialize(const std::vector<double> & vect, size_t & index) {
  fields_->deserialize(vect, index);
}
// -----------------------------------------------------------------------------
void IncrementQG::print(std::ostream & os) const {
  os << std::endl << "  Valid time: " << validTime();
  os << *fields_;
}
// -----------------------------------------------------------------------------
oops::LocalIncrement IncrementQG::getLocal(const GeometryQGIterator & iter) const {
  return fields_->getLocal(iter);
}
// -----------------------------------------------------------------------------
void IncrementQG::setLocal(const oops::LocalIncrement & values, const GeometryQGIterator & iter) {
  fields_->setLocal(values, iter);
}
// -----------------------------------------------------------------------------

}  // namespace qg
