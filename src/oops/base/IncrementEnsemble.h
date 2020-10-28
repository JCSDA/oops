/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_INCREMENTENSEMBLE_H_
#define OOPS_BASE_INCREMENTENSEMBLE_H_

#include <memory>
#include <string>
#include <utility>
#include <vector>

#include <boost/ptr_container/ptr_vector.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/base/Accumulator.h"
#include "oops/base/LinearVariableChangeBase.h"
#include "oops/base/StateEnsemble.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// \brief Ensemble of inrements
template<typename MODEL> class IncrementEnsemble {
  typedef Geometry<MODEL>            Geometry_;
  typedef Increment<MODEL>           Increment_;
  typedef LinearVariableChangeBase<MODEL>  LinearVariableChangeBase_;
  typedef State<MODEL>               State_;
  typedef StateEnsemble<MODEL>       StateEnsemble_;

  typedef typename boost::ptr_vector<LinearVariableChangeBase_> ChvarVec_;
  typedef typename ChvarVec_::const_reverse_iterator ircst_;

 public:
  /// Constructor
  IncrementEnsemble(const Geometry_ & resol, const Variables & vars,
                    const util::DateTime &, const int rank);
  /// \brief construct ensemble of perturbations as \p ens - \p mean; holding
  //         \p vars variables
  IncrementEnsemble(const StateEnsemble_ & ens, const State_ & mean, const Variables & vars);
  IncrementEnsemble(const eckit::Configuration &, const State_ &, const State_ &,
                    const Geometry_ &, const Variables &);

  /// Accessors
  size_t size() const {return ensemblePerturbs_.size();}
  Increment_ & operator[](const int ii) {return ensemblePerturbs_[ii];}
  const Increment_ & operator[](const int ii) const {return ensemblePerturbs_[ii];}

  /// Control variables
  const Variables & controlVariables() const {return vars_;}

  /// Release / reset
  void releaseMember();
  void appendMember(const Increment_ &);

 private:
  const Variables vars_;
  std::vector<Increment_> ensemblePerturbs_;
};

// ====================================================================================

template<typename MODEL>
IncrementEnsemble<MODEL>::IncrementEnsemble(const Geometry_ & resol, const Variables & vars,
                                            const util::DateTime & tslot, const int rank)
  : vars_(vars), ensemblePerturbs_()
{
  ensemblePerturbs_.reserve(rank);
  for (int m = 0; m < rank; ++m) {
    ensemblePerturbs_.emplace_back(resol, vars_, tslot);
  }
  Log::trace() << "IncrementEnsemble:contructor done" << std::endl;
}

// ====================================================================================

template<typename MODEL>
IncrementEnsemble<MODEL>::IncrementEnsemble(const StateEnsemble_ & ensemble,
                                            const State_ & mean, const Variables & vars)
  : vars_(vars), ensemblePerturbs_()
{
  ensemblePerturbs_.reserve(ensemble.size());
  for (size_t ii = 0; ii < ensemble.size(); ++ii) {
    ensemblePerturbs_.emplace_back(ensemble[ii].geometry(), vars, ensemble[ii].validTime());
    ensemblePerturbs_[ii].diff(ensemble[ii], mean);
  }
  Log::trace() << "IncrementEnsemble:contructor(StateEnsemble) done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
IncrementEnsemble<MODEL>::IncrementEnsemble(const eckit::Configuration & conf,
                                            const State_ & xb, const State_ & fg,
                                            const Geometry_ & resol, const Variables & vars)
  : vars_(vars), ensemblePerturbs_()
{
  // Get rank from config
  std::vector<eckit::LocalConfiguration> memberConfig;
  conf.get("members", memberConfig);

  // Check sizes and fill in timeslots
  util::DateTime tslot = xb.validTime();

  // Read inflation field
  std::unique_ptr<Increment_> inflationField;
  if (conf.has("inflation field")) {
    const eckit::LocalConfiguration inflationConfig(conf, "inflation field");
    inflationField.reset(new Increment_(resol, vars, tslot));
    inflationField->read(inflationConfig);
  }

  // Get inflation value
  double inflationValue = conf.getDouble("inflation value", 1.0);

  // Setup change of variable
  ChvarVec_ chvars;
  std::vector<eckit::LocalConfiguration> chvarconfs;
  conf.get("variable changes", chvarconfs);
  for (const auto & conf : chvarconfs) {
    chvars.push_back(LinearVariableChangeFactory<MODEL>::create(xb, fg, resol, conf));
  }
  // TODO(Benjamin): one change of variable for each timeslot

  // Read ensemble
  StateEnsemble_ ensemble(resol, conf);
  State_ bgmean = ensemble.mean();

  ensemblePerturbs_.reserve(ensemble.size());
  for (unsigned int ie = 0; ie < ensemble.size(); ++ie) {
    // Ensemble will be centered around ensemble mean
    Increment_ dx(resol, vars_, tslot);
    dx.diff(ensemble[ie], bgmean);

    // Apply inflation
    if (conf.has("inflation field")) {
      dx.schur_product_with(*inflationField);
    }
    dx *= inflationValue;

    // Apply inverse of the linear balance operator
    // K_1^{-1} K_2^{-1} .. K_N^{-1}
    for (ircst_ it = chvars.rbegin(); it != chvars.rend(); ++it) {
      dx = it->multiplyInverse(dx);
    }

    ensemblePerturbs_.emplace_back(std::move(dx));
  }
  Log::trace() << "IncrementEnsemble:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IncrementEnsemble<MODEL>::releaseMember() {
  ensemblePerturbs_.erase(ensemblePerturbs_.begin());
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IncrementEnsemble<MODEL>::appendMember(const Increment_ & dx) {
  ensemblePerturbs_.emplace_back(dx);
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_INCREMENTENSEMBLE_H_
