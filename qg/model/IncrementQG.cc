/*
 * (C) Copyright 2009-2016 ECMWF.
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

#include "eckit/config/LocalConfiguration.h"
#include "model/GomQG.h"
#include "model/LocationsQG.h"
#include "model/ModelBiasIncrement.h"
#include "model/ErrorCovarianceQG.h"
#include "model/FieldsQG.h"
#include "model/GeometryQG.h"
#include "model/StateQG.h"
#include "model/VariablesQG.h"
#include "oops/generic/UnstructuredGrid.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/Logger.h"

using oops::Log;

namespace qg {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
IncrementQG::IncrementQG(const GeometryQG & resol, const VariablesQG & vars,
                         const util::DateTime & vt)
  : fields_(new FieldsQG(resol, vars, vt)), stash_()
{
  fields_->zero();
  Log::trace() << "IncrementQG constructed." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementQG::IncrementQG(const GeometryQG & resol, const IncrementQG & other)
  : fields_(new FieldsQG(*other.fields_, resol)), stash_()
{
  Log::trace() << "IncrementQG constructed from other." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementQG::IncrementQG(const IncrementQG & other, const bool copy)
  : fields_(new FieldsQG(*other.fields_, copy)), stash_()
{
  Log::trace() << "IncrementQG copy-created." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementQG::IncrementQG(const IncrementQG & other)
  : fields_(new FieldsQG(*other.fields_)), stash_()
{
  Log::trace() << "IncrementQG copy-created." << std::endl;
}
// -----------------------------------------------------------------------------
IncrementQG::~IncrementQG() {
  Log::trace() << "IncrementQG destructed" << std::endl;
}
// -----------------------------------------------------------------------------
void IncrementQG::activateModel() {
// Should get variables from model. YT
  eckit::LocalConfiguration modelvars;
  modelvars.set("variables", "tl");
  VariablesQG vars(modelvars);
// Should get variables from model. YT
  stash_.reset(new FieldsQG(*fields_, vars));
  swap(fields_, stash_);
  ASSERT(fields_);
  ASSERT(stash_);
  Log::trace() << "IncrementQG activated for TLM" << std::endl;
}
// -----------------------------------------------------------------------------
void IncrementQG::deactivateModel() {
  swap(fields_, stash_);
  *fields_ = *stash_;
  stash_.reset();
  ASSERT(fields_);
  ASSERT(!stash_);
  Log::trace() << "IncrementQG deactivated for TLM" << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
void IncrementQG::diff(const StateQG & x1, const StateQG & x2) {
  ASSERT(this->validTime() == x1.validTime());
  ASSERT(this->validTime() == x2.validTime());
  Log::debug() << "IncrementQG:diff incr " << *fields_ << std::endl;
  Log::debug() << "IncrementQG:diff x1 " << x1.fields() << std::endl;
  Log::debug() << "IncrementQG:diff x2 " << x2.fields() << std::endl;
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
/// Interpolate to observation location
// -----------------------------------------------------------------------------
void IncrementQG::interpolateTL(const LocationsQG & locs, GomQG & cols) const {
  Log::debug() << "IncrementQG::interpolateTL fields in" << *fields_ << std::endl;
  fields_->interpolateTL(locs, cols);
  Log::debug() << "IncrementQG::interpolateTL fields out" << *fields_ << std::endl;
  Log::debug() << "IncrementQG::interpolateTL gom " << cols << std::endl;
}
// -----------------------------------------------------------------------------
void IncrementQG::interpolateAD(const LocationsQG & locs, const GomQG & cols) {
  Log::debug() << "IncrementQG::interpolateAD gom " << cols << std::endl;
  Log::debug() << "IncrementQG::interpolateAD fields in" << *fields_ << std::endl;
  fields_->interpolateAD(locs, cols);
  Log::debug() << "IncrementQG::interpolateAD fields out" << *fields_ << std::endl;
}
// -----------------------------------------------------------------------------
/// Convert to/from unstructured grid
// -----------------------------------------------------------------------------
void IncrementQG::convert_to(oops::UnstructuredGrid & ug) const {
  fields_->convert_to(ug);
}
// -----------------------------------------------------------------------------
void IncrementQG::convert_from(const oops::UnstructuredGrid & ug) {
  fields_->convert_from(ug);
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
void IncrementQG::print(std::ostream & os) const {
  os << std::endl << "  Valid time: " << validTime();
  os << *fields_;
}
// -----------------------------------------------------------------------------

}  // namespace qg
