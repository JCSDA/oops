/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_STATEENSEMBLE_H_
#define OOPS_BASE_STATEENSEMBLE_H_

#include <string>
#include <vector>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "oops/assimilation/Increment4D.h"
#include "oops/assimilation/State4D.h"
#include "oops/base/Accumulator.h"
#include "oops/base/LinearVariableChangeBase.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/util/DateTime.h"
#include "oops/util/dot_product.h"
#include "oops/util/Duration.h"
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Ensemble

template<typename MODEL> class StateEnsemble {
  typedef LinearVariableChangeBase<MODEL>  LinearVariableChangeBase_;
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;
  typedef State4D<MODEL>             State4D_;
  typedef Increment<MODEL>           Increment_;
  typedef Increment4D<MODEL>         Increment4D_;

  typedef typename boost::ptr_vector<LinearVariableChangeBase_> ChvarVec_;
  typedef typename ChvarVec_::const_reverse_iterator ircst_;

 public:
/// Constructor
  StateEnsemble();
  StateEnsemble(const Geometry_ & resol,
                const Variables & vars,
                const std::vector<util::DateTime> &,
                const int rank);
  StateEnsemble(const std::vector<util::DateTime> &,
                const eckit::Configuration &);

/// Destructor
  virtual ~StateEnsemble() {}

  /// Accessors
  unsigned int size() const {
    return rank_;
  }
  Increment4D_ & operator[](const int ii) {
    return ensemblePerturbs_[ii];
  }
  const Increment4D_ & operator[](const int ii) const {
    return ensemblePerturbs_[ii];
  }

  void linearize(const State4D_ &, const State4D_ &, const Geometry_ &);

  const Variables & controlVariables() const {return vars_;}

 private:
  const eckit::LocalConfiguration config_;

  unsigned int rank_;
  const std::vector<util::DateTime> timeslots_;
  const Variables vars_;
  boost::scoped_ptr<const Geometry_> resol_;

  boost::ptr_vector<Increment4D_> ensemblePerturbs_;
};

// ====================================================================================

template<typename MODEL>
StateEnsemble<MODEL>::StateEnsemble()
  : rank_(0), resol_(), ensemblePerturbs_()
{
  Log::trace() << "StateEnsemble:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
StateEnsemble<MODEL>::StateEnsemble(const Geometry_ & resol,
                                    const Variables & vars,
                                    const std::vector<util::DateTime> & timeslots,
                                    const int rank)
  : rank_(rank), timeslots_(timeslots),
    vars_(vars), resol_(), ensemblePerturbs_()
{
  resol_.reset(new Geometry_(resol));
  for (unsigned m = 0; m < rank_; ++m) {
    Increment4D_ * dx = new Increment4D_(*resol_, vars_, timeslots_);
    ensemblePerturbs_.push_back(dx);
  }
  Log::trace() << "StateEnsemble:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
StateEnsemble<MODEL>::StateEnsemble(const std::vector<util::DateTime> & timeslots,
                                    const eckit::Configuration & conf)
  : config_(conf), rank_(0), timeslots_(timeslots),
    vars_(conf), resol_(), ensemblePerturbs_()
{
  // Get rank from config
  std::vector<eckit::LocalConfiguration> memberConfig;
  config_.get("members", memberConfig);
  rank_ = memberConfig.size();
  Log::trace() << "StateEnsemble:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void StateEnsemble<MODEL>::linearize(const State4D_ & xb, const State4D_ & fg,
                                     const Geometry_ & resol) {
  // Check sizes
  ASSERT(xb.size() == timeslots_.size());
  ASSERT(fg.size() == timeslots_.size());
  for (unsigned jsub = 0; jsub < timeslots_.size(); ++jsub) {
     ASSERT(xb[jsub].validTime() == timeslots_[jsub]);
     ASSERT(fg[jsub].validTime() == timeslots_[jsub]);
  }
  std::vector<eckit::LocalConfiguration> memberConfig;
  config_.get("members", memberConfig);
  ASSERT(memberConfig.size() == rank_);

  // Set resolution
  resol_.reset(new Geometry_(resol));

  // Setup change of variable
  ChvarVec_ chvars;
  if (config_.has("variable_changes")) {
    std::vector<eckit::LocalConfiguration> chvarconfs;
    config_.get("variable_changes", chvarconfs);
    for (const auto & conf : chvarconfs) {
      chvars.push_back(LinearVariableChangeFactory<MODEL>::create(xb[0], fg[0], resol, conf));
    }
  }
  // TODO(Benjamin): one change of variable for each timeslot

  // Read ensemble
  std::vector<State4D_> ensemble;
  for (unsigned int ie = 0; ie < rank_; ++ie) {
    boost::scoped_ptr<State4D_> xx;
    if (memberConfig[ie].has("state")) {
      xx.reset(new State4D_(memberConfig[ie], vars_, *resol_));
    } else {
      State_ xx3D(*resol_, vars_, memberConfig[ie]);
      xx.reset(new State4D_(xx3D));
    }
    ensemble.push_back((*xx));
  }

  // Compute ensemble mean
  Accumulator<MODEL, State4D_, State4D_> bgmean(ensemble[0]);
  for (unsigned int ie = 0; ie < rank_; ++ie) {
    const double rr = 1.0/static_cast<double>(rank_);
    bgmean.accumul(rr, ensemble[ie]);
  }

  const double rk = 1.0 / sqrt((static_cast<double>(rank_) - 1.0));
  for (unsigned int ie = 0; ie < rank_; ++ie) {
    // Ensemble will be centered around ensemble mean
    Increment4D_ dx(*resol_, vars_, timeslots_);
    dx.diff(ensemble[ie], bgmean);

    // Apply inverse of the linear balance operator
    for (unsigned jsub = 0; jsub < timeslots_.size(); ++jsub) {
      // K_1^{-1} K_2^{-1} .. K_N^{-1}
      for (ircst_ it = chvars.rbegin(); it != chvars.rend(); ++it) {
        Increment_ dxchvarout = it->multiplyInverse(dx[jsub]);
        dx[jsub] = dxchvarout;
      }
    }

    Increment4D_ * dxunbalptr = new Increment4D_(dx);
    ensemblePerturbs_.push_back(dxunbalptr);

    // Rescale
    ensemblePerturbs_[ie] *= rk;
  }
}

// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_STATEENSEMBLE_H_
