/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef LORENZ95_OBSTABLE_H_
#define LORENZ95_OBSTABLE_H_

#include <fstream>
#include <map>
#include <ostream>
#include <string>
#include <vector>

#include "oops/interface/ObsSpaceBase.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"

namespace eckit {
  class Configuration;
}

namespace lorenz95 {
  class LocsL95;
  class ObsVec1D;

/// A Simple Observation Data Handler
/*!
 *  ObsTable defines a simple observation handler
 *  that mimicks the interfaces required from ODB.
 */

// -----------------------------------------------------------------------------
class ObsTable : public oops::ObsSpaceBase,
                 private util::ObjectCounter<ObsTable> {
 public:
  static const std::string classname() {return "lorenz95::ObsTable";}

  ObsTable(const eckit::Configuration &, const util::DateTime &, const util::DateTime &);
  ~ObsTable();

  void putdb(const std::string &, const std::vector<double> &) const;
  void getdb(const std::string &, std::vector<double> &) const;

  LocsL95 * locations(const util::DateTime & t1, const util::DateTime & t2) const;
  std::vector<int> timeSelect(const util::DateTime &, const util::DateTime &) const;
  void generateDistribution(const eckit::Configuration &);
  void printJo(const ObsVec1D &, const ObsVec1D &);
  unsigned int nobs() const {return times_.size();}

 private:
  void print(std::ostream &) const;
  void otOpen(const std::string &);
  void otWrite(const std::string &) const;

  const util::DateTime winbgn_;
  const util::DateTime winend_;

  std::vector<util::DateTime> times_;
  std::vector<double> locations_;
  mutable std::map<std::string, std::vector<double> > data_;

  std::string nameIn_;
  std::string nameOut_;
};
// -----------------------------------------------------------------------------
}  // namespace lorenz95

#endif  // LORENZ95_OBSTABLE_H_
