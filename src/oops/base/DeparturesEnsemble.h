/*
 * (C) Copyright 2018 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_BASE_DEPARTURESENSEMBLE_H_
#define OOPS_BASE_DEPARTURESENSEMBLE_H_

#include <memory>
#include <vector>

#include "oops/base/Departures.h"
#include "oops/base/ObsEnsemble.h"
#include "oops/base/Observations.h"
#include "oops/base/ObsSpaces.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Ensemble

template<typename MODEL> class DeparturesEnsemble {
  typedef Departures<MODEL>          Departures_;
  typedef ObsEnsemble<MODEL>         ObsEnsemble_;
  typedef Observations<MODEL>        Observations_;
  typedef ObsSpaces<MODEL>           ObsSpaces_;

 public:
/// Constructor
  explicit DeparturesEnsemble(const ObsEnsemble_ &);
  DeparturesEnsemble(const ObsEnsemble_ &, const Observations_ &);
  DeparturesEnsemble(const ObsSpaces_ &, const DeparturesEnsemble &);

/// Destructor
  virtual ~DeparturesEnsemble() {}

  /// Accessors
  unsigned int size() const {
    return ensemblePerturbs_.size();
  }
  Departures_ & operator[](const int ii) {
    return *ensemblePerturbs_[ii];
  }
  const Departures_ & operator[](const int ii) const {
    return *ensemblePerturbs_[ii];
  }

 private:
  std::vector<std::shared_ptr<Departures_>> ensemblePerturbs_;
};

// ====================================================================================

template<typename MODEL>
DeparturesEnsemble<MODEL>::DeparturesEnsemble(const ObsEnsemble_ & ens)
  : DeparturesEnsemble(ens, ens.mean()) {}

// -----------------------------------------------------------------------------

template<typename MODEL>
DeparturesEnsemble<MODEL>::DeparturesEnsemble(const ObsEnsemble_ & ens,
                                              const Observations_ & ensmean)
  : ensemblePerturbs_()
{
  for (unsigned i = 0; i < ens.size(); ++i) {
    std::shared_ptr<Departures_> y(new Departures_(ens[i] - ensmean));
    ensemblePerturbs_.push_back(y);
  }
  Log::trace() << "DeparturesEnsemble:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
DeparturesEnsemble<MODEL>::DeparturesEnsemble(const ObsSpaces_ & obsdb,
                                              const DeparturesEnsemble & other)
  : ensemblePerturbs_() {
  for (std::size_t jj = 0; jj < other.size(); ++jj) {
    std::shared_ptr<Departures_> tmp(new Departures_(obsdb, other[jj]));
    ensemblePerturbs_.push_back(tmp);
  }
  Log::trace() << "Local DeparturesEnsemble created" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_DEPARTURESENSEMBLE_H_
