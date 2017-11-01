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

#include "util/Logger.h"
#include "eckit/config/Configuration.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/abor1_cpp.h"

#include "lorenz95/FieldL95.h"
#include "lorenz95/GomL95.h"
#include "lorenz95/IncrementL95.h"
#include "lorenz95/LocsL95.h"
#include "lorenz95/ModelBias.h"
#include "lorenz95/ModelL95.h"
#include "lorenz95/ModelTrajectory.h"
#include "lorenz95/Resolution.h"



using oops::Log;

namespace lorenz95 {

// -----------------------------------------------------------------------------
/// Constructor, destructor
// -----------------------------------------------------------------------------
StateL95::StateL95(const Resolution & resol, const NoVariables &,
                   const util::DateTime & vt)
  : fld_(resol), time_(vt)
{
  Log::trace() << "StateL95::StateL95 created" << std::endl;
}
// -----------------------------------------------------------------------------
StateL95::StateL95(const Resolution & resol, const eckit::Configuration & file)
  : fld_(resol), time_(util::DateTime())
{
  this->read(file);
  Log::trace() << "StateL95::StateL95 created and read in." << std::endl;
}
// -----------------------------------------------------------------------------
StateL95::StateL95(const Resolution & resol, const StateL95 & xx)
  : fld_(resol), time_(xx.time_)
{
  fld_ = xx.fld_;
  Log::trace() << "StateL95::StateL95 created by interpolation." << std::endl;
}
// -----------------------------------------------------------------------------
StateL95::StateL95(const StateL95 & xx)
  : fld_(xx.fld_), time_(xx.time_)
{
  Log::trace() << "StateL95::StateL95 copy-created." << std::endl;
}
// -----------------------------------------------------------------------------
StateL95::~StateL95() {
  Log::trace() << "StateL95::StateL95 destructed." << std::endl;
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
/// Interpolate to observation location
// -----------------------------------------------------------------------------
void StateL95::interpolate(const LocsL95 & locs, GomL95 & vals) const {
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
  Log::trace() << "StateL95::read opening " << filename << std::endl;
  std::ifstream fin(filename.c_str());
  if (!fin.is_open()) ABORT("StateL95::read: Error opening file");

  int resol;
  fin >> resol;
  ASSERT(fld_.resol() == resol);

  std::string stime;
  fin >> stime;
  const util::DateTime tt(stime);
  const util::DateTime tc(config.getString("date"));
  if (tc != tt) {
    ABORT("StateL95::read: date and data file inconsistent.");
  }
  time_ = tt;

  fld_.read(fin);

  fin.close();
  Log::trace() << "StateL95::read: file closed." << std::endl;
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

  Log::trace() << "StateL95::write opening " << filename << std::endl;
  std::ofstream fout(filename.c_str());
  if (!fout.is_open()) ABORT("StateL95::write: Error opening file");

  fout << fld_.resol() << std::endl;
  fout << time_ << std::endl;
  fld_.write(fout);
  fout << std::endl;

  fout.close();
  Log::trace() << "StateL95::write file closed." << std::endl;
}
// -----------------------------------------------------------------------------
void StateL95::print(std::ostream & os) const {
  os << std::endl << " Valid time: " << time_;
  os << std::endl << fld_;
}
// -----------------------------------------------------------------------------

}  // namespace lorenz95
