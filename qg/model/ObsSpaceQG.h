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

#include "eckit/geometry/Point2.h"
#include "eckit/mpi/Comm.h"

#include "oops/base/ObsSpaceBase.h"
#include "oops/base/Variables.h"
#include "oops/util/DateTime.h"

#include "oops/qg/LocationsQG.h"
#include "oops/qg/QgFortran.h"

namespace eckit {
  class Configuration;
}

namespace qg {
  class ObsVecQG;

/// \brief ObsSpace for QG model
//  \details ObsSpaceQG is created for each obs type. The underlying Fortran
//  structure (key_) is created for each matching input-output filename pair
//  (i.e. different obstypes can be stored in the same Fortran structure).
//  For mapping between ObsSpaceQG and Fortran structures,
//  ObsSpaceQG::theObsFileRegister_ map is used
class ObsSpaceQG : public oops::ObsSpaceBase {
 public:
  /// create full ObsSpace (read or generate data)
  ObsSpaceQG(const eckit::Configuration &, const eckit::mpi::Comm &,
             const util::DateTime &, const util::DateTime &);
  /// create local ObsSpace
  ObsSpaceQG(const ObsSpaceQG &, const eckit::geometry::Point2 &,
             const eckit::Configuration &);
  ~ObsSpaceQG();

  /// read data or metadata
  void getdb(const std::string &, int &) const;
  /// save data or metadata
  void putdb(const std::string &, const int &) const;

  /// check if variable is in ObsSpace
  bool has(const std::string & col) const;

  /// create locations between times (\p t1, \p t2]
  std::unique_ptr<LocationsQG> locations(const util::DateTime & t1,
                               const util::DateTime & t2) const;

  void printJo(const ObsVecQG &, const ObsVecQG &) const;

  /// return number of observations (unique locations)
  int nobs() const;

  /// return variables simulated by ObsOperators
  const oops::Variables & obsvariables() const { return obsvars_; }

  /// observation type
  const std::string & obsname() const {return obsname_;}

  /// local observations indices
  const std::vector<int> & localobs() const { return localobs_;}

  /// interface with Fortran
  const F90odb & toFortran() const {return key_;}

 private:
  void print(std::ostream &) const;

  F90odb key_;                       // pointer to Fortran structure
  const std::string obsname_;        // corresponds with obstype
  const util::DateTime winbgn_;      // window for the observations
  const util::DateTime winend_;
  oops::Variables obsvars_;          // variables simulated by ObsOperators
  std::vector<int> localobs_;        // indices of local observations
  bool isLocal_;                     // true if it's a local subset
  const eckit::mpi::Comm & comm_;    // MPI communicator associated with ObsSpace

  // defines mapping for Fortran structures
  static std::map < std::string, F90odb > theObsFileRegister_;
  static int theObsFileCount_;  // number of files used
};

}  // namespace qg

#endif  // QG_MODEL_OBSSPACEQG_H_
