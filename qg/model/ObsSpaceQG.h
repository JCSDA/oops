/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef QG_MODEL_OBSSPACEQG_H_
#define QG_MODEL_OBSSPACEQG_H_

#include <map>
#include <memory>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/mpi/Comm.h"

#include "oops/base/ObsSpaceBase.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"
#include "oops/util/parameters/OptionalParameter.h"
#include "oops/util/parameters/Parameters.h"
#include "oops/util/parameters/RequiredParameter.h"

#include "oops/qg/LocationsQG.h"
#include "oops/qg/ObsIteratorQG.h"
#include "oops/qg/QgFortran.h"

namespace qg {
  class ObsIteratorQG;

// -----------------------------------------------------------------------------
/// Contents of the `engine` YAML section.
class ObsDataParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsDataParameters, Parameters)

 public:
  /// File path and file type
  oops::RequiredParameter<std::string> obsfile{"obsfile", this};
};

/// Contents of the `obsdatain` or `obsdataout` YAML section.
class ObsEngineParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsEngineParameters, Parameters)

 public:
  /// File path.
  oops::RequiredParameter<ObsDataParameters> engine{"engine", this};
};

/// Options specifying locations of artificial observations.
class ObsLocGenParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsLocGenParameters, Parameters)

 public:
  oops::RequiredParameter<std::vector<float>> lon{"lon",
       "longitudes for obs locations, degrees between [-180.0, 180.0]", this};
  oops::RequiredParameter<std::vector<float>> lat{"lat",
       "latitudes for obs locations, degrees between [5.0, 84.0]", this};
  oops::RequiredParameter<std::vector<float>> z{"z",
       "depths for obs locations, meters between 0 and 1.e4", this};
};

/// Options controlling generation of artificial observations.
class ObsGenerateParameters : public oops::Parameters {
  OOPS_CONCRETE_PARAMETERS(ObsGenerateParameters, Parameters)

 public:
  oops::RequiredParameter<util::Duration> begin{"begin", this};
  oops::RequiredParameter<util::Duration> obsPeriod{"obs period", this};
  /// Number of observations to generate in each time slot.
  oops::OptionalParameter<int> obsDensity{"obs density", this};
  /// Or locations of observations to generate in each time slot.
  oops::OptionalParameter<ObsLocGenParameters> obsLocs{"obs locations", this};

  oops::RequiredParameter<int> nval{"nval", this};
  oops::RequiredParameter<double> obsError{"obs error", this};
};

/// \brief Configuration parameters for the QG model's ObsSpace.
class ObsSpaceQGParameters : public oops::ObsSpaceParametersBase {
  OOPS_CONCRETE_PARAMETERS(ObsSpaceQGParameters, ObsSpaceParametersBase)

 public:
  /// Type of observations.
  oops::RequiredParameter<std::string> obsType{"obs type", this};
  /// File from which to load observations.
  oops::OptionalParameter<ObsEngineParameters> obsdatain{"obsdatain", this};
  /// File to which to save observations and analysis.
  oops::OptionalParameter<ObsEngineParameters> obsdataout{"obsdataout", this};
  /// Options controlling generation of artificial observations.
  oops::OptionalParameter<ObsGenerateParameters> generate{"generate", this};
};

// -----------------------------------------------------------------------------
/// \brief ObsSpace for QG model
/// \details ObsSpaceQG is created for each obs type. The underlying Fortran
/// structure (key_) is created for each matching input-output filename pair
/// (i.e. different obstypes can be stored in the same Fortran structure).
/// For mapping between ObsSpaceQG and Fortran structures,
/// ObsSpaceQG::theObsFileRegister_ map is used
class ObsSpaceQG : public oops::ObsSpaceBase {
 public:
  typedef ObsSpaceQGParameters Parameters_;

  /// create full ObsSpace (read or generate data)
  ObsSpaceQG(const Parameters_ &, const eckit::mpi::Comm &,
             const util::DateTime &, const util::DateTime &, const eckit::mpi::Comm &);
  ~ObsSpaceQG();

  /// save and close file
  void save() const;

  /// read data or metadata
  void getdb(const std::string &, int &) const;
  /// save data or metadata
  void putdb(const std::string &, const int &) const;

  /// sample the location of each observation in the time window with a single path
  std::unique_ptr<LocationsQG> locations() const;

  /// return number of observations (unique locations)
  int nobs() const;

  /// return variables to be processed
  const oops::Variables & obsvariables() const { return obsvars_; }

  /// return variables simulated by ObsOperators
  const oops::Variables & assimvariables() const { return assimvars_; }

  /// observation type
  const std::string & obsname() const {return obsname_;}

  /// iterator to the first observation
  ObsIteratorQG begin() const;
  /// iterator to the observation past-the-last
  ObsIteratorQG end() const;

  /// interface with Fortran
  const F90odb & toFortran() const {return key_;}

 private:
  void print(std::ostream &) const;

  mutable F90odb key_;               // pointer to Fortran structure
  const std::string obsname_;        // corresponds with obstype
  const util::DateTime winbgn_;      // window for the observations
  const util::DateTime winend_;
  oops::Variables assimvars_;          // variables simulated by ObsOperators
  oops::Variables obsvars_;          // variables that are observed

  // defines mapping for Fortran structures
  static std::map < std::string, F90odb > theObsFileRegister_;
  static int theObsFileCount_;  // number of files used
};

}  // namespace qg

#endif  // QG_MODEL_OBSSPACEQG_H_
