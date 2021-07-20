/*
 * (C) Copyright 2018-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_DEPARTURESENSEMBLE_H_
#define OOPS_BASE_DEPARTURESENSEMBLE_H_

#include <Eigen/Dense>
#include <memory>
#include <vector>

#include "oops/base/Departures.h"
#include "oops/base/ObsSpaces.h"
#include "oops/interface/ObsDataVector.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Ensemble of Departures (can hold ensemble perturbations in the observation space)
template<typename OBS> class DeparturesEnsemble {
  typedef Departures<OBS>          Departures_;
  typedef ObsSpaces<OBS>           ObsSpaces_;
  template <typename DATA> using ObsData_ = ObsDataVector<OBS, DATA>;
  template <typename DATA> using ObsDataVec_ = std::vector<std::shared_ptr<ObsData_<DATA>>>;

 public:
  /// Creates ensemble of empty Departures size \p nens
  DeparturesEnsemble(const ObsSpaces_ &, const size_t nens);

  /// Accessors and size
  size_t size() const {return ensemblePerturbs_.size();}
  Departures_ & operator[](const size_t ii) {return ensemblePerturbs_[ii];}
  const Departures_ & operator[](const size_t ii) const {return ensemblePerturbs_[ii];}

/// pack ensemble of dep. as contiguous block of memory
  Eigen::MatrixXd packEigen(const Departures_ &) const;

 private:
  std::vector<Departures_> ensemblePerturbs_;   // ensemble perturbations
};

// ====================================================================================

template<typename OBS>
DeparturesEnsemble<OBS>::DeparturesEnsemble(const ObsSpaces_ & obsdb, const size_t nens)
      : ensemblePerturbs_() {
  ensemblePerturbs_.reserve(nens);
  for (size_t iens = 0; iens < nens; ++iens) {
    ensemblePerturbs_.emplace_back(obsdb);
  }
  Log::trace() << "DeparturesEnsemble created" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename OBS>
Eigen::MatrixXd DeparturesEnsemble<OBS>::packEigen(const Departures_ & mask) const {
  std::size_t myNobs = ensemblePerturbs_[0].packEigenSize(mask);
  std::size_t myNens = ensemblePerturbs_.size();

  Eigen::MatrixXd depEns(myNens, myNobs);
  for (std::size_t iens = 0; iens < myNens; ++iens) {
    depEns.row(iens) = ensemblePerturbs_[iens].packEigen(mask);
  }
  Log::trace() << "DeparturesEnsemble::packEigen() completed" << std::endl;
  return depEns;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_DEPARTURESENSEMBLE_H_
