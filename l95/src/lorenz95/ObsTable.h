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
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

namespace lorenz95 {
  class ObsIterator;

// -----------------------------------------------------------------------------
/// Contents of the `engine` YAML section.
class ObsDataParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsDataParameters, Parameters)

 public:
  /// File path and file type
  oops::RequiredParameter<std::string> obsfile{"obsfile", this};
};
// -----------------------------------------------------------------------------
/// Contents of the `obsdatain` or `obsdataout` YAML section.
class ObsEngineParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsEngineParameters, Parameters)

 public:
  /// File path.
  oops::RequiredParameter<ObsDataParameters> engine{"engine", this};
};
// -----------------------------------------------------------------------------
/// Options controlling generation of artificial observations.
class ObsGenerateParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsGenerateParameters, Parameters)

 public:
  oops::RequiredParameter<util::Duration> obsFrequency{"obs_frequency", this};
  /// Number of observations to generate in each time slot.
  oops::RequiredParameter<int> obsDensity{"obs_density", this};
  oops::RequiredParameter<double> obsError{"obs_error", this};
};

// -----------------------------------------------------------------------------
/// \brief Configuration parameters for the L95 model's ObsSpace.
class ObsTableParameters : public oops::ObsSpaceParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsTableParameters, ObsSpaceParametersBase)

 public:
  /// File from which to load observations.
  oops::OptionalParameter<ObsEngineParameters> obsdatain{"obsdatain", this};
  /// File to which to save observations and analysis.
  oops::OptionalParameter<ObsEngineParameters> obsdataout{"obsdataout", this};
  /// Options controlling generation of artificial observations.
  oops::OptionalParameter<ObsGenerateParameters> generate{"generate", this};
};

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

  typedef ObsTableParameters Parameters_;

  ObsTable(const Parameters_ &, const eckit::mpi::Comm &,
           const util::DateTime &, const util::DateTime &, const eckit::mpi::Comm &);
  ~ObsTable();

  void save() const;

  void putdb(const std::string &, const std::vector<int> &) const;
  void putdb(const std::string &, const std::vector<float> &) const;
  void putdb(const std::string &, const std::vector<double> &) const;
  void getdb(const std::string &, std::vector<int> &) const;
  void getdb(const std::string &, std::vector<float> &) const;
  void getdb(const std::string &, std::vector<double> &) const;

  bool has(const std::string & col) const;
  void generateDistribution(const ObsGenerateParameters & params);
  void random(std::vector<double> &) const;
  unsigned int nobs() const {return times_.size();}
  const std::vector<double> & locations() const { return locations_; }
  const std::vector<util::DateTime> & times() const { return times_; }
  const oops::Variables & obsvariables() const { return obsvars_; }
  const oops::Variables & assimvariables() const { return assimvars_; }
  const std::string & obsname() const {return obsname_;}

  /// iterator to the first observation
  ObsIterator begin() const;
  /// iterator to the observation past-the-last
  ObsIterator end() const;

 private:
  void print(std::ostream &) const;
  void otOpen(const std::string &);
  void otWrite(const std::string &) const;

  const util::DateTime winbgn_;
  const util::DateTime winend_;

  std::vector<util::DateTime> times_;
  std::vector<double> locations_;
  mutable std::map<std::string, std::vector<double> > data_;

  const eckit::mpi::Comm & comm_;
  const oops::Variables obsvars_;
  const oops::Variables assimvars_;
  std::string nameIn_;
  std::string nameOut_;
  const std::string obsname_ = "Lorenz 95";
};
// -----------------------------------------------------------------------------
}  // namespace lorenz95

#endif  // LORENZ95_OBSTABLE_H_
