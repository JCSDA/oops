/*
 * (C) Copyright 2018-2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_DEPARTURESENSEMBLE_H_
#define OOPS_BASE_DEPARTURESENSEMBLE_H_

#include <vector>

#include "oops/base/Departures.h"
#include "oops/base/ObsSpaces.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Ensemble of Departures (can hold ensemble perturbations in the observation space)
template<typename MODEL> class DeparturesEnsemble {
  typedef Departures<MODEL>          Departures_;
  typedef ObsSpaces<MODEL>           ObsSpaces_;
 public:
  /// Creates ensemble of empty Departures size \p nens
  DeparturesEnsemble(const ObsSpaces_ &, const size_t nens);
  /// Creates ensemble of local Departures from full Departures \p other based on local
  /// observations \p local
  DeparturesEnsemble(const ObsSpaces_ & local, const DeparturesEnsemble & other);

  /// Accessors and size
  size_t size() const {return ensemblePerturbs_.size();}
  Departures_ & operator[](const size_t ii) {return ensemblePerturbs_[ii];}
  const Departures_ & operator[](const size_t ii) const {return ensemblePerturbs_[ii];}

 private:
  std::vector<Departures_> ensemblePerturbs_;   // ensemble perturbations
};

// ====================================================================================

template<typename MODEL>
DeparturesEnsemble<MODEL>::DeparturesEnsemble(const ObsSpaces_ & obsdb, const size_t nens)
      : ensemblePerturbs_() {
  ensemblePerturbs_.reserve(nens);
  for (size_t iens = 0; iens < nens; ++iens) {
    ensemblePerturbs_.emplace_back(obsdb);
  }
  Log::trace() << "DeparturesEnsemble created" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
DeparturesEnsemble<MODEL>::DeparturesEnsemble(const ObsSpaces_ & local,
                                              const DeparturesEnsemble & other)
      : ensemblePerturbs_() {
  ensemblePerturbs_.reserve(other.size());
  for (const auto & dep : other.ensemblePerturbs_) {
    ensemblePerturbs_.emplace_back(local, dep);
  }
  Log::trace() << "Local DeparturesEnsemble created" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_DEPARTURESENSEMBLE_H_
