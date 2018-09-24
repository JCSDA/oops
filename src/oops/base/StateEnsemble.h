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
  typedef Increment<MODEL>           Increment_;

  typedef typename boost::ptr_vector<LinearVariableChangeBase_> ChvarVec_;
  typedef typename ChvarVec_::const_reverse_iterator ircst_;

 public:
/// Constructor
  StateEnsemble(const util::DateTime &, const eckit::Configuration &);

/// Destructor
  virtual ~StateEnsemble() {}

  /// Accessors
  unsigned int size() const {
    return rank_;
  }
  Increment_ & operator[](const int ii) {
    return ensemblePerturbs_[ii];
  }
  const Increment_ & operator[](const int ii) const {
    return ensemblePerturbs_[ii];
  }

  void linearize(const State_ &, const State_ &, const Geometry_ &);

  const Variables & controlVariables() const {return vars_;}

 private:
  const eckit::LocalConfiguration config_;

  unsigned int rank_;
  const util::DateTime validTime_;
  const Variables vars_;
  boost::scoped_ptr<const Geometry_> resol_;

  boost::ptr_vector<Increment_> ensemblePerturbs_;
};

// ====================================================================================

template<typename MODEL>
StateEnsemble<MODEL>::StateEnsemble(const util::DateTime & validTime,
                                    const eckit::Configuration & conf)
  : config_(conf), rank_(conf.getInt("members")), validTime_(validTime),
    vars_(eckit::LocalConfiguration(conf, "variables")),
    resol_(), ensemblePerturbs_()
{
  Log::trace() << "StateEnsemble:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void StateEnsemble<MODEL>::linearize(const State_ & xb, const State_ & fg,
                                     const Geometry_ & resol) {
  ASSERT(xb.validTime() == validTime_);
  resol_.reset(new Geometry_(resol));

// Setup change of variable
  ChvarVec_ chvars;
  if (config_.has("variable_changes")) {
    std::vector<eckit::LocalConfiguration> chvarconfs;
    config_.get("variable_changes", chvarconfs);
    for (const auto & conf : chvarconfs) {
      chvars.push_back(LinearVariableChangeFactory<MODEL>::create(xb, fg, resol, conf));
    }
  }

// Read ensemble and compute mean
  State_ xblr(*resol_, xb);
  Accumulator<MODEL, State_, State_> bgmean(*resol_, vars_, validTime_);
  const double rr = 1.0/static_cast<double>(rank_);

  std::vector<eckit::LocalConfiguration> confs;
  config_.get("state", confs);
  ASSERT(confs.size() == rank_);

  State_ xread(xblr);
  std::vector<State_> ensemble;
  for (unsigned int jm = 0; jm < rank_; ++jm) {
    xread.read(confs[jm]);
    ASSERT(xread.validTime() == validTime_);
    ensemble.push_back(xread);

//  Compute ensemble mean
    bgmean.accumul(rr, xread);
  }

  const double rk = 1.0 / sqrt((static_cast<double>(rank_) - 1.0));
  for (unsigned int jm = 0; jm < rank_; ++jm) {
//  Ensemble will be centered around ensemble mean
    boost::scoped_ptr<Increment_> dx(new Increment_(*resol_, vars_, validTime_));
    dx->diff(ensemble[jm], bgmean);

//  Apply inverse of the linear balance operator

    // K_1^{-1} K_2^{-1} .. K_N^{-1}
    for (ircst_ it = chvars.rbegin(); it != chvars.rend(); ++it) {
      Increment_ dxchvarout = it->multiplyInverse(*dx);
      dx.reset(new Increment_(dxchvarout));
    }

    Increment_ * dxunbalptr = new Increment_(*dx);
    ensemblePerturbs_.push_back(dxunbalptr);

//  Rescale
    ensemblePerturbs_[jm] *= rk;
  }
}

// -----------------------------------------------------------------------------
}  // namespace oops

#endif  // OOPS_BASE_STATEENSEMBLE_H_
