/*
 * (C) Copyright 2018-2025 UCAR
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
  size_t size() const {return nens_;}

  /// pack ensemble of dep. as contiguous block of memory for a given localization
  Eigen::MatrixXf packEigen(const Departures_ &);

  /// getData: get data of imem-th member in ensemblePerturbs_ and store in a Departures
  /// setData: set data of imem-th member in ensemblePerturbs_ from a Departures
  /// Note: To modify data of imem-th member in ensemblePerturbs_,
  ///       e.g., ensemblePerturbs_.getData(imem) += dep,
  ///             or ensemblePerturbs_.getData(imem).mask(mask_),
  ///       it will
  ///       first declare a temporary Departures,
  ///       call getData to get data of imem-th member and store in a Deartures,
  ///       apply operations on the Departures,
  ///       and call setData to store updated data in the Departuers back to
  ///       imem-th member in ensemblePerturbs_.
  ///       In the example of operator += above for imem-th member,
  ///       the code is like
  ///           Departures_ tmpDeps(this->obspaces_);
  ///           tmpDeps = ensemblePerturbs_.getData(imem);
  ///           tmpDeps += dep;
  ///           ensemblePerturbs_.setData(imem, tmpDeps);
  void setData(const size_t imem, const Departures_ &);
  Departures_ getData(const size_t imem) const;

 private:
  Eigen::MatrixXf ensemblePerturbs_;
  const ObsSpaces_ & obsdb_;
  std::size_t nens_;
  /// A vector of serial sizes at each observation space
  std::vector<size_t> serialSizes_;
};
// ====================================================================================
template<typename OBS>
DeparturesEnsemble<OBS>::DeparturesEnsemble(const ObsSpaces_ & obsdb, const size_t nens)
  : ensemblePerturbs_(), obsdb_(obsdb), nens_(nens)
{
  Departures_ dep(obsdb_);
  serialSizes_ = dep.serialSizes();
  const size_t ntot = std::accumulate(serialSizes_.begin(), serialSizes_.end(), 0);
  ensemblePerturbs_.resize(nens, ntot);
  Log::trace() << "DeparturesEnsemble created" << std::endl;
}

// -----------------------------------------------------------------------------
template<typename OBS>
Eigen::MatrixXf DeparturesEnsemble<OBS>::packEigen(const Departures_ & mask) {
  std::vector<std::vector<size_t>> indices = this->getData(0).maskAndSerialIndices(mask);

  size_t myNobs = 0;
  for (size_t ii = 0; ii != indices.size(); ++ii) {
      myNobs += indices[ii].size();
  }

  const size_t myNens = nens_;
  Eigen::MatrixXf depEns(myNens, myNobs);

  size_t istart = 0;
  size_t icol = 0;
  for (size_t ii = 0; ii != serialSizes_.size(); ++ii) {
    for (size_t jj = 0; jj != indices[ii].size(); ++jj) {
        depEns.col(icol) = ensemblePerturbs_.col(indices[ii][jj] + istart);
        icol++;
    }
    istart += serialSizes_[ii];
  }

  Log::trace() << "DeparturesEnsemble::packEigen() completed" << std::endl;
  return depEns;
}
// -----------------------------------------------------------------------------
template<typename OBS>
void DeparturesEnsemble<OBS>::setData(const size_t imem,
                              const Departures_ & dep) {
  std::vector<double> valvec;
  for (size_t ii = 0; ii != serialSizes_.size(); ++ii) {
    dep[ii].serialize(valvec);
  }
  const Eigen::VectorXd valvec_eigen = Eigen::Map<Eigen::VectorXd>(valvec.data(),
                                                                   valvec.size());
  this->ensemblePerturbs_.row(imem) = valvec_eigen.cast<float>();
}
// -----------------------------------------------------------------------------
template<typename OBS>
Departures<OBS> DeparturesEnsemble<OBS>::getData(const size_t imem) const {
  Departures_ dep(obsdb_);
  const Eigen::VectorXd valvec_eigen = ensemblePerturbs_.row(imem).template cast<double>();
  std::vector<double> valvec(valvec_eigen.data(),
                             valvec_eigen.data() + valvec_eigen.size());
  size_t istart = 0;
  for (size_t ii = 0; ii != serialSizes_.size(); ++ii) {
    dep[ii].deserialize(valvec, istart);
  }
  return dep;
}
// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_DEPARTURESENSEMBLE_H_
