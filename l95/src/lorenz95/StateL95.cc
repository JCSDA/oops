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
  const std::string &exp = parameters.exp;
  const std::string &type = parameters.type;
  std::string filename = dir+"/"+exp+"."+type;

  if (type == "ens") {
    if (parameters.member.value() == boost::none)
      throw eckit::BadValue("'member' was not set in the parameters passed to write() "
                            "even though 'type' was set to '" + type + "'", Here());
    const int &memb = *parameters.member.value();
    filename += "."+std::to_string(memb);
  }

  if (type == "fc" || type == "ens") {
    if (parameters.date.value() == boost::none)
      throw eckit::BadValue("'date' was not set in the parameters passed to write() "
                            "even though 'type' was set to '" + type + "'", Here());
    const util::DateTime &antime = *parameters.date.value();
    filename += "."+antime.toString();
    const util::Duration step = time_ - antime;
    filename += "."+step.toString();
  }

  if (type == "an") {
    filename += "."+time_.toString();
  }

  if (type == "krylov") {
    if (parameters.iteration.value() == boost::none)
      throw eckit::BadValue("'iteration' was not set in the parameters passed to write() "
                            "even though 'type' was set to '" + type + "'", Here());
    const int &iter = *parameters.iteration.value();
    filename += "."+std::to_string(iter)+"."+time_.toString();
  }

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
  size_t nn = 3;
  nn += fld_.serialSize();
  nn += time_.serialSize();
  return nn;
}
// -----------------------------------------------------------------------------
void StateL95::serialize(std::vector<double> & vect) const {
  vect.push_back(1001.0);
  fld_.serialize(vect);
  vect.push_back(2002.0);
  time_.serialize(vect);
  vect.push_back(3003.0);
}
// -----------------------------------------------------------------------------
void StateL95::deserialize(const std::vector<double> & vect, size_t & index) {
  size_t ii = index + this->serialSize();
  ASSERT(vect.at(index) == 1001.0);
  ++index;
  fld_.deserialize(vect, index);
  ASSERT(vect.at(index) == 2002.0);
  ++index;
  time_.deserialize(vect, index);
  ASSERT(vect.at(index) == 3003.0);
  ++index;
  ASSERT(index == ii);
}
// -----------------------------------------------------------------------------

}  // namespace lorenz95
