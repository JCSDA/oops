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
#include <vector>

#include <boost/ptr_container/ptr_vector.hpp>

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
#include "oops/util/Logger.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Ensemble
template<typename MODEL> class IncrementEnsemble {
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
  IncrementEnsemble(const Geometry_ & resol,
                    const Variables & vars,
                    const std::vector<util::DateTime> &,
                    const int rank);
  IncrementEnsemble(const eckit::Configuration &,
                    const State4D_ &, const State4D_ &, const Geometry_ &);

  /// Destructor
  virtual ~IncrementEnsemble() {}

  /// Accessors
  unsigned int size() const {
    return rank_;
  }
  Increment4D_ & operator[](const int ii) {
    return *ensemblePerturbs_[ii];
  }
  const Increment4D_ & operator[](const int ii) const {
    return *ensemblePerturbs_[ii];
  }

  /// Control variables
  const Variables & controlVariables() const {return vars_;}

  /// Release / reset
  void releaseMember();
  void resetMember(const Increment4D_ &);

 private:
  unsigned int rank_;
  const Variables vars_;
  std::vector<std::unique_ptr<Increment4D_>> ensemblePerturbs_;
};

// ====================================================================================

template<typename MODEL>
IncrementEnsemble<MODEL>::IncrementEnsemble(const Geometry_ & resol,
                                            const Variables & vars,
                                            const std::vector<util::DateTime> & timeslots,
                                            const int rank)
  : rank_(rank), vars_(vars), ensemblePerturbs_()
{
  for (unsigned m = 0; m < rank_; ++m) {
    Increment4D_ * dx = new Increment4D_(resol, vars_, timeslots);
    ensemblePerturbs_.push_back(std::unique_ptr<Increment4D_>(dx));
  }
  Log::trace() << "IncrementEnsemble:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
IncrementEnsemble<MODEL>::IncrementEnsemble(const eckit::Configuration & conf,
                                            const State4D_ & xb, const State4D_ & fg,
                                            const Geometry_ & resol)
  : rank_(0), vars_(conf), ensemblePerturbs_()
{
  // Get rank from config
  std::vector<eckit::LocalConfiguration> memberConfig;
  conf.get("members", memberConfig);
  rank_ = memberConfig.size();

  // Check sizes and fill in timeslots
  ASSERT(xb.size() == fg.size());
  std::vector<util::DateTime> timeslots(xb.size());
  for (unsigned jsub = 0; jsub < xb.size(); ++jsub) {
     ASSERT(xb[jsub].validTime() == fg[jsub].validTime());
     timeslots[jsub] = xb[jsub].validTime();
  }

  // Setup change of variable
  ChvarVec_ chvars;
  if (conf.has("variable_changes")) {
    std::vector<eckit::LocalConfiguration> chvarconfs;
    conf.get("variable_changes", chvarconfs);
    for (const auto & conf : chvarconfs) {
      chvars.push_back(LinearVariableChangeFactory<MODEL>::create(xb[0], fg[0], resol, conf));
    }
  }
  // TODO(Benjamin): one change of variable for each timeslot

  // Read ensemble
  std::vector<State4D_> ensemble;
  for (unsigned int ie = 0; ie < rank_; ++ie) {
    std::unique_ptr<State4D_> xx;
    if (memberConfig[ie].has("state")) {
      xx.reset(new State4D_(memberConfig[ie], vars_, resol));
    } else {
      State_ xx3D(resol, vars_, memberConfig[ie]);
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
    Increment4D_ dx(resol, vars_, timeslots);
    dx.diff(ensemble[ie], bgmean);

    // Apply inverse of the linear balance operator
    for (unsigned jsub = 0; jsub < timeslots.size(); ++jsub) {
      // K_1^{-1} K_2^{-1} .. K_N^{-1}
      for (ircst_ it = chvars.rbegin(); it != chvars.rend(); ++it) {
        Increment_ dxchvarout = it->multiplyInverse(dx[jsub]);
        dx[jsub] = dxchvarout;
      }
    }

    Increment4D_ * dxunbalptr = new Increment4D_(dx);
    ensemblePerturbs_.push_back(std::unique_ptr<Increment4D_>(dxunbalptr));

    // Rescale
    *ensemblePerturbs_[ie] *= rk;
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
void IncrementEnsemble<MODEL>::resetMember(const Increment4D_ & dx) {
  Increment4D_ * dxptr = new Increment4D_(dx);
  ensemblePerturbs_.push_back(std::unique_ptr<Increment4D_>(dxptr));
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_INCREMENTENSEMBLE_H_
