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

#include "eckit/mpi/Comm.h"

#include "oops/base/ObsSpaceBase.h"
#include "oops/base/ObsVariables.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"

namespace eckit {
  class Configuration;
}

namespace lorenz95 {
  class ObsIterator;

// -----------------------------------------------------------------------------
/// A Simple Observation Data Handler
/*!
 *  ObsTable defines a simple observation handler
 *  that mimicks the interfaces required from ODB.
 */
class ObsTable : public oops::ObsSpaceBase,
                 private util::ObjectCounter<ObsTable> {
 public:
  static const std::string classname() {return "lorenz95::ObsTable";}

  ObsTable(const eckit::Configuration &, const eckit::mpi::Comm &,
           const util::TimeWindow &, const eckit::mpi::Comm &);
  ~ObsTable();

  void save() const;
  void append(const std::string & appendDir);

  void putdb(const std::string &, const std::vector<int> &) const;
  void putdb(const std::string &, const std::vector<float> &) const;
  void putdb(const std::string &, const std::vector<double> &) const;
  void getdb(const std::string &, std::vector<int> &) const;
  void getdb(const std::string &, std::vector<float> &) const;
  void getdb(const std::string &, std::vector<double> &) const;

  bool has(const std::string & col) const;
  void generateDistribution(const eckit::Configuration &);
  void random(std::vector<double> &) const;
  unsigned int nobs() const {return times_.size();}
  const std::vector<double> & locations() const { return locations_; }
  const std::vector<util::DateTime> & times() const { return times_; }
  const oops::ObsVariables & obsvariables() const { return obsvars_; }
  const oops::ObsVariables & assimvariables() const { return assimvars_; }
  const std::string & obsname() const {return obsname_;}

  /// iterator to the first observation
  ObsIterator begin() const;
  /// iterator to the observation past-the-last
  ObsIterator end() const;

 private:
  void print(std::ostream &) const;
  void otOpen(const std::string &);
  void otWrite(const std::string &) const;

  std::vector<util::DateTime> times_;
  std::vector<double> locations_;
  mutable std::map<std::string, std::vector<double> > data_;

  const eckit::mpi::Comm & comm_;
  const util::TimeWindow timeWindow_;
  const oops::ObsVariables obsvars_;
  const oops::ObsVariables assimvars_;
  std::string nameIn_;
  std::string nameOut_;
  const std::string obsname_ = "Lorenz 95";
};
// -----------------------------------------------------------------------------
}  // namespace lorenz95

#endif  // LORENZ95_OBSTABLE_H_
