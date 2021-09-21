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
#include "oops/base/Geometry.h"
#include "oops/base/Increment.h"
#include "oops/base/LinearVariableChangeBase.h"
#include "oops/base/State.h"
#include "oops/base/StateEnsemble.h"
#include "oops/base/Variables.h"
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
  IncrementEnsemble(const eckit::Configuration &, const State_ &, const State_ &,
                    const Geometry_ &, const Variables &);
  /// \brief construct ensemble of perturbations by reading them from disk
  IncrementEnsemble(const Geometry_ &, const Variables &, const eckit::Configuration &);
  /// \brief construct ensemble of perturbations by reading two state ensembles (one member at a
  //         time) and taking the  difference of each set of pairs
  IncrementEnsemble(const Geometry_ &, const Variables &, const eckit::Configuration &,
                    const eckit::Configuration &);

  void write(const eckit::Configuration &) const;

  /// Accessors
  size_t size() const {return ensemblePerturbs_.size();}
  Increment_ & operator[](const int ii) {return ensemblePerturbs_[ii];}
  const Increment_ & operator[](const int ii) const {return ensemblePerturbs_[ii];}

  /// Control variables
  const Variables & controlVariables() const {return vars_;}

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

// -----------------------------------------------------------------------------

template<typename MODEL>
IncrementEnsemble<MODEL>::IncrementEnsemble(const eckit::Configuration & conf,
                                            const State_ & xb, const State_ & fg,
                                            const Geometry_ & resol, const Variables & vars)
  : vars_(vars), ensemblePerturbs_()
{
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
IncrementEnsemble<MODEL>::IncrementEnsemble(const Geometry_ & resol, const Variables & vars,
                                            const eckit::Configuration & config)
  : vars_(vars), ensemblePerturbs_()
{
  std::vector<eckit::LocalConfiguration> memberConfig;
  config.get("members", memberConfig);

  // Datetime for ensemble
  util::DateTime tslot = util::DateTime(config.getString("date"));

  // Reserve memory to hold ensemble
  ensemblePerturbs_.reserve(memberConfig.size());

  // Loop over all ensemble members
  for (size_t jj = 0; jj < memberConfig.size(); ++jj) {
    Increment_ dx(resol, vars_, tslot);
    dx.read(memberConfig[jj]);
    ensemblePerturbs_.emplace_back(std::move(dx));
  }
  Log::trace() << "IncrementEnsemble:contructor (by reading increment ensemble) done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
IncrementEnsemble<MODEL>::IncrementEnsemble(const Geometry_ & resol, const Variables & vars,
                                            const eckit::Configuration & configBase,
                                            const eckit::Configuration & configPert)
  : vars_(vars), ensemblePerturbs_()
{
  std::vector<eckit::LocalConfiguration> memberConfigBase;
  configBase.get("members", memberConfigBase);

  std::vector<eckit::LocalConfiguration> memberConfigPert;
  configPert.get("members", memberConfigPert);

  // Ensure input ensembles are of the same size
  ASSERT(memberConfigBase.size() == memberConfigPert.size());

  // Reserve memory to hold ensemble
  ensemblePerturbs_.reserve(memberConfigBase.size());

  // Loop over all ensemble members
  for (size_t jj = 0; jj < memberConfigBase.size(); ++jj) {
    State_ xBase(resol, memberConfigBase[jj]);
    State_ xPert(resol, memberConfigPert[jj]);
    Increment_ dx(resol, vars_, xBase.validTime());
    dx.diff(xBase, xPert);
    ensemblePerturbs_.emplace_back(std::move(dx));
  }
  Log::trace() << "IncrementEnsemble:contructor (by diffing state ensembles) done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void IncrementEnsemble<MODEL>::write(const eckit::Configuration & config) const
{
  eckit::LocalConfiguration outConfig(config);
  for (size_t ii=0; ii < size(); ++ii) {
    outConfig.set("member", ii+1);
    ensemblePerturbs_[ii].write(outConfig);
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_INCREMENTENSEMBLE_H_
