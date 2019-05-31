/*
 * (C) Copyright 2009-2016 ECMWF.
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
#include "lorenz95/FieldL95.h"
#include "lorenz95/GomL95.h"
#include "lorenz95/IncrementL95.h"
#include "lorenz95/LocsL95.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelL95.h"
#include "lorenz95/ModelTrajectory.h"
#include "lorenz95/Nothing.h"
#include "lorenz95/Resolution.h"
#include "oops/generic/UnstructuredGrid.h"
#include "oops/util/abor1_cpp.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {
  class Variables;
}

namespace lorenz95 {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
StateL95::StateL95(const Resolution & resol, const oops::Variables &,
                   const util::DateTime & vt)
  : fld_(resol), time_(vt)
{
  oops::Log::trace() << "StateL95::StateL95 created" << std::endl;
}
// -----------------------------------------------------------------------------
StateL95::StateL95(const Resolution & resol, const oops::Variables &,
                   const eckit::Configuration & conf)
  : fld_(resol), time_(conf.getString("date"))
{
  oops::Log::trace() << "StateL95::StateL95 conf " << conf << std::endl;
  if (conf.has("filename")) {
    this->read(conf);
  } else {
    fld_.generate(conf);
  }
  oops::Log::trace() << "StateL95::StateL95 created and read in." << std::endl;
}
// -----------------------------------------------------------------------------
StateL95::StateL95(const Resolution & resol, const StateL95 & xx)
  : fld_(resol), time_(xx.time_)
{
  fld_ = xx.fld_;
  oops::Log::trace() << "StateL95::StateL95 created by interpolation." << std::endl;
}
// -----------------------------------------------------------------------------
StateL95::StateL95(const StateL95 & xx)
  : fld_(xx.fld_), time_(xx.time_)
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
  return *this;
}
// -----------------------------------------------------------------------------
/// Get state values at obs locations
// -----------------------------------------------------------------------------
void StateL95::getValues(const LocsL95 & locs, const oops::Variables &, GomL95 & vals) const {
  fld_.interp(locs, vals);
}
// -----------------------------------------------------------------------------
void StateL95::getValues(const LocsL95 & locs, const oops::Variables &, GomL95 & vals,
                           Nothing &) const {
  fld_.interp(locs, vals);
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
  const std::string filename(config.getString("filename"));
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
  std::string dir = config.getString("datadir");
  std::string exp = config.getString("exp");
  std::string type = config.getString("type");
  std::string filename = dir+"/"+exp+"."+type;

  if (type == "ens") {
    std::string memb = config.getString("member");
    filename += "."+memb;
  }

  if (type == "fc" || type == "ens") {
    const util::DateTime antime(config.getString("date"));
    filename += "."+antime.toString();
    const util::Duration step = time_ - antime;
    filename += "."+step.toString();
  }

  if (type == "an") {
    filename += "."+time_.toString();
  }

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
oops::GridPoint StateL95::getPoint(const Iterator & i) const {
  std::vector<std::string> vars;
  vars.push_back("x");
  std::vector<double> vals;
  vals.push_back(fld_[i.index()]);
  std::vector<int> varlens;
  varlens.push_back(1);
  return oops::GridPoint(oops::Variables(vars), vals, varlens);
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


}  // namespace lorenz95
