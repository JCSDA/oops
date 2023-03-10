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

#include "lorenz95/StateL95.h"

#include <fstream>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "eckit/exception/Exceptions.h"

#include "lorenz95/FieldL95.h"
#include "lorenz95/GomL95.h"
#include "lorenz95/IncrementL95.h"
#include "lorenz95/LocsL95.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelL95.h"
#include "lorenz95/ModelTrajectory.h"
#include "lorenz95/Resolution.h"

#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"
#include "oops/util/stringFunctions.h"

namespace oops {
  class Variables;
}
namespace sf = util::stringfunctions;

namespace lorenz95 {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
StateL95::StateL95(const Resolution & resol, const oops::Variables & vars,
                   const util::DateTime & vt)
  : fld_(resol), time_(vt), vars_(vars)
{
  oops::Log::trace() << "StateL95::StateL95 created" << std::endl;
}
// -----------------------------------------------------------------------------
StateL95::StateL95(const Resolution & resol, const Parameters_ & parameters)
  : fld_(resol), time_(parameters.date), vars_({"x"})
{
  oops::Log::trace() << "StateL95::StateL95 conf " << parameters << std::endl;
  if (parameters.filename.value() != boost::none) {
    this->read(parameters);
  } else if (parameters.analyticInit.value() != boost::none) {
    fld_.generate(*parameters.analyticInit.value());
  }
  oops::Log::trace() << "StateL95::StateL95 created and read in." << std::endl;
}
// -----------------------------------------------------------------------------
StateL95::StateL95(const Resolution & resol, const StateL95 & xx)
  : fld_(resol), time_(xx.time_), vars_(xx.vars_)
{
  fld_ = xx.fld_;
  oops::Log::trace() << "StateL95::StateL95 created by interpolation." << std::endl;
}
// -----------------------------------------------------------------------------
StateL95::StateL95(const StateL95 & xx)
  : fld_(xx.fld_), time_(xx.time_), vars_(xx.vars_)
{
  oops::Log::trace() << "StateL95::StateL95 copy-created." << std::endl;
}
// -----------------------------------------------------------------------------
StateL95::~StateL95() {
  oops::Log::trace() << "StateL95::StateL95 destructed." << std::endl;
}
// -----------------------------------------------------------------------------
/// Basic operators
// -----------------------------------------------------------------------------
StateL95 & StateL95::operator=(const StateL95 & rhs) {
  fld_ = rhs.fld_;
  time_ = rhs.time_;
  vars_ = rhs.vars_;
  return *this;
}
// -----------------------------------------------------------------------------
/// Interactions with Increments
// -----------------------------------------------------------------------------
StateL95 & StateL95::operator+=(const IncrementL95 & dx) {
  ASSERT(time_ == dx.validTime());
  fld_ += dx.getField();
  return *this;
}
// -----------------------------------------------------------------------------
/// Atlas
// -----------------------------------------------------------------------------
void StateL95::toFieldSet(atlas::FieldSet & fset) const {
  ASSERT(fset.empty());

  // TODO(FH): revisit how the Field is created when updating to "version 2" of atlas interfaces:
  // The "proper" way to create the Field would be via the FunctionSpace's createField, so that the
  // Field is linked to the FunctionSpace. This may be needed for atlas's haloExchange method.
  const int resol = fld_.resol();
  atlas::Field fld("field", atlas::array::DataType::real64(), atlas::array::ArrayShape({resol, 1}));

  auto view = atlas::array::make_view<double, 2>(fld);
  for (size_t jj = 0; jj < static_cast<size_t>(resol); ++jj) {
    view(jj, 0) = fld_[jj];
  }
  fset.add(fld);
}
// -----------------------------------------------------------------------------
void StateL95::fromFieldSet(const atlas::FieldSet & fset) {
  ASSERT(!fset.empty());

  const int resol = fld_.resol();
  const auto & view = atlas::array::make_view<double, 2>(fset["field"]);
  for (size_t jj = 0; jj < static_cast<size_t>(resol); ++jj) {
    fld_[jj] = view(jj, 0);
  }
}
// -----------------------------------------------------------------------------
/// Utilities
// -----------------------------------------------------------------------------
void StateL95::read(const Parameters_ & parameters) {
  std::string filename(parameters.filename.value().value());
  sf::swapNameMember(parameters.member.value(), filename);
  oops::Log::trace() << "StateL95::read opening " << filename << std::endl;
  std::ifstream fin(filename.c_str());
  if (!fin.is_open()) ABORT("StateL95::read: Error opening file: " + filename);

  int resol;
  fin >> resol;
  ASSERT(fld_.resol() == resol);

  std::string stime;
  fin >> stime;
  const util::DateTime tt(stime);
  if (time_ != tt) {
    ABORT("StateL95::read: date and data file inconsistent.");
  }

  fld_.read(fin);

  fin.close();
  oops::Log::trace() << "StateL95::read: file closed." << std::endl;
}
// -----------------------------------------------------------------------------
void StateL95::write(const WriteParameters_ & parameters) const {
  const std::string &dir = parameters.datadir;
  const std::string &type = *parameters.type.value();
  std::string filename = dir+"/"+*parameters.prefix.value();

  if (type == "krylov") {
    if (parameters.iteration.value() == boost::none)
      throw eckit::BadValue("'iteration' was not set in the parameters passed to write() "
                            "even though 'type' was set to '" + type + "'", Here());
    const int &iter = *parameters.iteration.value();
    filename += "."+std::to_string(iter)+"."+time_.toString();
  }

  filename += ".l95";

  sf::swapNameMember(parameters.member.value(), filename);

  oops::Log::trace() << "StateL95::write opening " << filename << std::endl;
  std::ofstream fout(filename.c_str());
  if (!fout.is_open()) ABORT("StateL95::write: Error opening file: " + filename);

  fout << fld_.resol() << std::endl;
  fout << time_ << std::endl;
  fld_.write(fout);
  fout << std::endl;

  fout.close();
  oops::Log::trace() << "StateL95::write file closed." << std::endl;
}
// -----------------------------------------------------------------------------
void StateL95::print(std::ostream & os) const {
  os << std::endl << " Valid time: " << time_;
  os << std::endl << fld_;
}
// -----------------------------------------------------------------------------
/// For accumulator
// -----------------------------------------------------------------------------
void StateL95::zero() {
  fld_.zero();
}
// -----------------------------------------------------------------------------
void StateL95::accumul(const double & zz, const StateL95 & xx) {
  fld_.axpy(zz, xx.fld_);
}
// -----------------------------------------------------------------------------
/// Serialize - deserialize
// -----------------------------------------------------------------------------
size_t StateL95::serialSize() const {
  size_t nn = 0;
  nn += fld_.serialSize();
  nn += time_.serialSize();
  return nn;
}
// -----------------------------------------------------------------------------
void StateL95::serialize(std::vector<double> & vect) const {
  fld_.serialize(vect);
  time_.serialize(vect);
}
// -----------------------------------------------------------------------------
void StateL95::deserialize(const std::vector<double> & vect, size_t & index) {
  const size_t ii = index + this->serialSize();
  fld_.deserialize(vect, index);
  time_.deserialize(vect, index);
  ASSERT(index == ii);
}
// -----------------------------------------------------------------------------

}  // namespace lorenz95
