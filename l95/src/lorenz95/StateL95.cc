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

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "lorenz95/FieldL95.h"
#include "lorenz95/GomL95.h"
#include "lorenz95/IncrementL95.h"
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
StateL95::StateL95(const Resolution & resol, const eckit::Configuration & conf)
  : fld_(resol), time_(conf.getString("date")), vars_({"x"})
{
  oops::Log::trace() << "StateL95::StateL95 conf " << conf << std::endl;
  if (conf.has("filename")) {
    this->read(conf);
  } else {
    fld_.generate(eckit::LocalConfiguration(conf, "analytic init"));
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
StateL95::StateL95(const oops::Variables & vars, const StateL95 & xx) : StateL95(xx) {}
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
void StateL95::read(const eckit::Configuration & config) {
  std::string filename(config.getString("filename"));
  sf::swapNameMember(config, filename);
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
void StateL95::write(const eckit::Configuration & config) const {
  const std::string dir = config.getString("datadir");
  const std::string type = config.getString("type");
  std::string filename = dir+"/"+config.getString("prefix");

  if (type == "krylov") {
    if (!config.has("iteration"))
      throw eckit::BadValue("'iteration' was not set in the parameters passed to write() "
                            "even though 'type' was set to '" + type + "'", Here());
    const int iter = config.getInt("iteration");
    filename += "."+std::to_string(iter)+"."+time_.toString();
  }

  filename += ".l95";
  sf::swapNameMember(config, filename);

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
