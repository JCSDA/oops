/*
 * (C) Copyright 2020 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#ifndef OOPS_ASSIMILATION_LETKFSOLVERBASE_H_
#define OOPS_ASSIMILATION_LETKFSOLVERBASE_H_

#include <map>
#include <memory>
#include <string>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Departures.h"
#include "oops/base/DeparturesEnsemble.h"
#include "oops/base/IncrementEnsemble.h"
#include "oops/base/ObsErrors.h"
#include "oops/interface/GeometryIterator.h"
#include "oops/util/Logger.h"

namespace oops {

/// \brief Base class for LETKF-type solvers
template <typename MODEL>
class LETKFSolverBase {
  typedef Departures<MODEL>         Departures_;
  typedef DeparturesEnsemble<MODEL> DeparturesEnsemble_;
  typedef GeometryIterator<MODEL>   GeometryIterator_;
  typedef IncrementEnsemble<MODEL>  IncrementEnsemble_;
  typedef ObsErrors<MODEL>          ObsErrors_;

 public:
  LETKFSolverBase() {}
  virtual ~LETKFSolverBase() = default;

  /// KF update + posterior inflation at a grid point location (GeometryIterator_)
  virtual void measurementUpdate(const Departures_ &, const DeparturesEnsemble_ &,
                          const ObsErrors_ &, const IncrementEnsemble_ &,
                          const GeometryIterator_ &, IncrementEnsemble_ &);

 private:
  /// Computes weights
  virtual void computeWeights(const Departures_ &, const DeparturesEnsemble_ &,
                              const ObsErrors_ &) = 0;
  /// Applies weights and adds posterior inflation
  virtual void applyWeights(const IncrementEnsemble_ &, IncrementEnsemble_ &,
                            const GeometryIterator_ &) = 0;
};

// -----------------------------------------------------------------------------
template <typename MODEL>
void LETKFSolverBase<MODEL>::measurementUpdate(const Departures_ & dy,
                                               const DeparturesEnsemble_ & Yb,
                                               const ObsErrors_ & R,
                                               const IncrementEnsemble_ & bkg_pert,
                                               const GeometryIterator_ & i,
                                               IncrementEnsemble_ & ana_pert) {
  if (dy.nobs() == 0) {
    // no obs. so no need to update Wa_ and wa_
    // ana_pert[i]=bkg_pert[i]
    for (size_t itime=0; itime < bkg_pert[0].size(); ++itime) {
      for (size_t iens=0; iens < bkg_pert.size(); ++iens) {
        LocalIncrement gp = bkg_pert[iens][itime].getLocal(i);
        ana_pert[iens][itime].setLocal(gp, i);
      }
    }
  } else {
    // if obs are present do normal KF update
    computeWeights(dy, Yb, R);
    applyWeights(bkg_pert, ana_pert, i);
  }
}

// =============================================================================

/// \brief factory for LETKF solvers
template <typename MODEL>
class LETKFSolverFactory {
 public:
  static std::unique_ptr<LETKFSolverBase<MODEL>> create(const eckit::Configuration &,
                                                        size_t);
  virtual ~LETKFSolverFactory() = default;
 protected:
  explicit LETKFSolverFactory(const std::string &);
 private:
  virtual LETKFSolverBase<MODEL> * make(const eckit::Configuration &, size_t) = 0;
  static std::map < std::string, LETKFSolverFactory<MODEL> * > & getMakers() {
    static std::map < std::string, LETKFSolverFactory<MODEL> * > makers_;
    return makers_;
  }
};

// -----------------------------------------------------------------------------

template<class MODEL, class T>
class LETKFSolverMaker : public LETKFSolverFactory<MODEL> {
  virtual LETKFSolverBase<MODEL> * make(const eckit::Configuration & conf, size_t nens)
    { return new T(conf, nens); }
 public:
  explicit LETKFSolverMaker(const std::string & name) : LETKFSolverFactory<MODEL>(name) {}
};

// =============================================================================

template <typename MODEL>
LETKFSolverFactory<MODEL>::LETKFSolverFactory(const std::string & name) {
  if (getMakers().find(name) != getMakers().end()) {
    Log::error() << name << " already registered in LETKF solver factory." << std::endl;
    ABORT("Element already registered in LETKFSolverFactory.");
  }
  getMakers()[name] = this;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
std::unique_ptr<LETKFSolverBase<MODEL>>
LETKFSolverFactory<MODEL>::create(const eckit::Configuration & conf, size_t nens) {
  Log::trace() << "LETKFSolverBase<MODEL>::create starting" << std::endl;
  const std::string id = conf.getString("letkf.solver");
  typename std::map<std::string, LETKFSolverFactory<MODEL>*>::iterator
    jloc = getMakers().find(id);
  if (jloc == getMakers().end()) {
    Log::error() << id << " does not exist in LETKF solver factory." << std::endl;
    Log::error() << "LETKF solver Factory has " << getMakers().size() << " elements:" << std::endl;
    for (typename std::map<std::string, LETKFSolverFactory<MODEL>*>::const_iterator
         jj = getMakers().begin(); jj != getMakers().end(); ++jj) {
       Log::error() << "A " << jj->first << " LETKFSolver" << std::endl;
    }
    ABORT("Element does not exist in LETKFSolverFactory.");
  }
  std::unique_ptr<LETKFSolverBase<MODEL>> ptr(jloc->second->make(conf, nens));
  Log::trace() << "LETKFSolverBase<MODEL>::create done" << std::endl;
  return ptr;
}

// -----------------------------------------------------------------------------

}  // namespace oops
#endif  // OOPS_ASSIMILATION_LETKFSOLVERBASE_H_
