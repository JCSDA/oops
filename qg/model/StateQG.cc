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
#include <utility>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

#include "model/FieldsQG.h"
#include "model/GeometryQG.h"
#include "model/GomQG.h"
#include "model/IncrementQG.h"
#include "model/LocationsQG.h"
#include "model/ModelBias.h"
#include "model/ModelQG.h"


namespace qg {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
StateQG::StateQG(const GeometryQG & resol, const oops::Variables & vars,
                 const util::DateTime & vt)
  : fields_(new FieldsQG(resol, vars, lbc_, vt))
{
  oops::Log::trace() << "StateQG::StateQG created." << std::endl;
}
// -----------------------------------------------------------------------------
StateQG::StateQG(const GeometryQG & resol, const eckit::Configuration & file)
  : fields_()
{
  oops::Variables vars({"x"});
  if (file.has("state variables")) vars = oops::Variables(file, "state variables");
  oops::Log::trace() << "StateQG::StateQG variables: " << vars << std::endl;
  fields_.reset(new FieldsQG(resol, vars, 1, util::DateTime()));
  if (file.has("analytic init")) {
    fields_->analytic_init(file);
  } else if (file.has("read_from_file")) {
    const int read_from_file = file.getInt("read_from_file");
    if (read_from_file == 0) {
       fields_->analytic_init(file);
    } else if (read_from_file == 1) {
      fields_->read(file);
    }
  } else {
    fields_->read(file);
  }

  ASSERT(fields_);
  oops::Log::trace() << "StateQG::StateQG created and read in." << std::endl;
}
// -----------------------------------------------------------------------------
StateQG::StateQG(const GeometryQG & resol, const StateQG & other)
  : fields_(new FieldsQG(*other.fields_, resol))
{
  ASSERT(fields_);
  oops::Log::trace() << "StateQG::StateQG created by interpolation." << std::endl;
}
// -----------------------------------------------------------------------------
StateQG::StateQG(const oops::Variables & vars, const StateQG & other)
  : fields_(new FieldsQG(*other.fields_, vars))
{
  ASSERT(fields_);
  oops::Log::trace() << "StateQG::StateQG created by copying selected variables." << std::endl;
}
// -----------------------------------------------------------------------------
StateQG::StateQG(const StateQG & other)
  : fields_(new FieldsQG(*other.fields_))
{
  ASSERT(fields_);
  oops::Log::trace() << "StateQG::StateQG copied." << std::endl;
}
// -----------------------------------------------------------------------------
StateQG::~StateQG() {
  oops::Log::trace() << "StateQG::StateQG destructed." << std::endl;
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
/// Interpolate full fields
// -----------------------------------------------------------------------------
void StateQG::changeResolution(const StateQG & other) {
  fields_->changeResolution(*other.fields_);
  oops::Log::trace() << "StateQG interpolated" << std::endl;
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
/// Serialization
// -----------------------------------------------------------------------------
size_t StateQG::serialSize() const {
  size_t nn = fields_->serialSize();
  return nn;
}
// -----------------------------------------------------------------------------
void StateQG::serialize(std::vector<double> & vect) const {
  fields_->serialize(vect);
}
// -----------------------------------------------------------------------------
void StateQG::deserialize(const std::vector<double> & vect, size_t & index) {
  fields_->deserialize(vect, index);
}
// -----------------------------------------------------------------------------
void StateQG::print(std::ostream & os) const {
  os << std::endl << "  Valid time: " << validTime();
  os << *fields_;
}
// -----------------------------------------------------------------------------
/// ATLAS FieldSet
// -----------------------------------------------------------------------------
void StateQG::toFieldSet(atlas::FieldSet & fset) const {
  fields_->toFieldSet(fset);
}
void StateQG::fromFieldSet(const atlas::FieldSet & fset) {
  fields_->fromFieldSet(fset);
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
