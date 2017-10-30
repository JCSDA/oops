/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#include "model/StateQG.h"

#include <algorithm>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "model/GomQG.h"
#include "model/LocationsQG.h"
#include "model/ModelBias.h"
#include "model/FieldsQG.h"
#include "model/GeometryQG.h"
#include "model/IncrementQG.h"
#include "model/ModelQG.h"
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
StateQG::StateQG(const GeometryQG & resol, const VariablesQG & vars,
                 const util::DateTime & vt)
  : fields_(new FieldsQG(resol, vars, vt)), stash_()
{
  Log::trace() << "StateQG::StateQG created." << std::endl;
}
// -----------------------------------------------------------------------------
StateQG::StateQG(const GeometryQG & resol, const eckit::Configuration & file)
  : fields_(), stash_()
{
// Should get variables from file. YT
  eckit::LocalConfiguration modelvars;
  modelvars.set("variables", "cv");
  VariablesQG vars(modelvars);
// Should get variables from file. YT
  fields_.reset(new FieldsQG(resol, vars, util::DateTime()));
  fields_->read(file);

  ASSERT(fields_);
  Log::trace() << "StateQG::StateQG created and read in." << std::endl;
}
// -----------------------------------------------------------------------------
StateQG::StateQG(const GeometryQG & resol, const StateQG & other)
  : fields_(new FieldsQG(*other.fields_, resol)), stash_()
{
  ASSERT(fields_);
  Log::trace() << "StateQG::StateQG created by interpolation." << std::endl;
}
// -----------------------------------------------------------------------------
StateQG::StateQG(const StateQG & other)
  : fields_(new FieldsQG(*other.fields_)), stash_()
{
  ASSERT(fields_);
  Log::trace() << "StateQG::StateQG copied." << std::endl;
}
// -----------------------------------------------------------------------------
StateQG::~StateQG() {
  Log::trace() << "StateQG::StateQG destructed." << std::endl;
}
// -----------------------------------------------------------------------------
void StateQG::activateModel() {
// Should get variables from model. YT
  eckit::LocalConfiguration modelvars;
  modelvars.set("variables", "nl");
  VariablesQG vars(modelvars);
// Should get variables from model. YT
  stash_.reset(new FieldsQG(*fields_, vars));
  swap(fields_, stash_);
  ASSERT(fields_);
  ASSERT(stash_);
  Log::trace() << "StateQG activated for Model" << std::endl;
}
// -----------------------------------------------------------------------------
void StateQG::deactivateModel() {
  swap(fields_, stash_);
  *fields_ = *stash_;
  stash_.reset();
  ASSERT(fields_);
  ASSERT(!stash_);
  Log::trace() << "StateQG deactivated for Model" << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
StateQG & StateQG::operator=(const StateQG & rhs) {
  ASSERT(fields_);
  *fields_ = *rhs.fields_;
  return *this;
}
// -----------------------------------------------------------------------------
/// Interpolate to observation location
// -----------------------------------------------------------------------------
void StateQG::interpolate(const LocationsQG & locs, GomQG & cols) const {
  fields_->interpolate(locs, cols);
}
// -----------------------------------------------------------------------------
/// Interpolate full fields
// -----------------------------------------------------------------------------
void StateQG::changeResolution(const StateQG & other) {
  fields_->changeResolution(*other.fields_);
  Log::trace() << "StateQG interpolated" << std::endl;
}
// -----------------------------------------------------------------------------
/// Interactions with Increments
// -----------------------------------------------------------------------------
StateQG & StateQG::operator+=(const IncrementQG & dx) {
  ASSERT(this->validTime() == dx.validTime());
  ASSERT(fields_);
  fields_->add(dx.fields());
  return *this;
}
// -----------------------------------------------------------------------------
/// Convert to/from unstructured grid
// -----------------------------------------------------------------------------
void StateQG::convert_to(oops::UnstructuredGrid & ug) const {
  fields_->convert_to(ug);
}
// -----------------------------------------------------------------------------
void StateQG::convert_from(const oops::UnstructuredGrid & ug) {
  fields_->convert_from(ug);
}
// -----------------------------------------------------------------------------
/// I/O and diagnostics
// -----------------------------------------------------------------------------
void StateQG::read(const eckit::Configuration & files) {
  fields_->read(files);
}
// -----------------------------------------------------------------------------
void StateQG::write(const eckit::Configuration & files) const {
  fields_->write(files);
}
// -----------------------------------------------------------------------------
void StateQG::print(std::ostream & os) const {
  os << std::endl << "  Valid time: " << validTime();
  os << *fields_;
}
// -----------------------------------------------------------------------------
/// For accumulator
// -----------------------------------------------------------------------------
void StateQG::zero() {
  fields_->zero();
}
// -----------------------------------------------------------------------------
void StateQG::accumul(const double & zz, const StateQG & xx) {
  fields_->axpy(zz, *xx.fields_);
}
// -----------------------------------------------------------------------------

}  // namespace qg
