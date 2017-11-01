/*
 * (C) Copyright 2009-2016 ECMWF.
 * 
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0. 
 * In applying this licence, ECMWF does not waive the privileges and immunities 
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_ENSEMBLE_H_
#define OOPS_BASE_ENSEMBLE_H_

#include <string>
#include <vector>

#include <boost/ptr_container/ptr_vector.hpp>
#include <boost/scoped_ptr.hpp>

#include "eckit/config/LocalConfiguration.h"
#include "util/Logger.h"
#include "oops/base/Accumulator.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Increment.h"
#include "oops/interface/State.h"
#include "oops/interface/Variables.h"
#include "util/DateTime.h"
#include "util/Duration.h"
#include "util/abor1_cpp.h"
#include "util/dot_product.h"

namespace oops {

// -----------------------------------------------------------------------------

/// Ensemble

template<typename MODEL> class Ensemble {
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;
  typedef Increment<MODEL>           Increment_;
  typedef Variables<MODEL>           Variables_;

 public:
/// Constructor
  Ensemble(const util::DateTime&, const eckit::Configuration&);

/// Destructor
  virtual ~Ensemble() {}

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

  void linearize(const State_ &, const Geometry_ &);

  const Variables_ & controlVariables() const {return vars_;}

 private:
  const eckit::LocalConfiguration config_;
  const util::DateTime validTime_;
  boost::scoped_ptr<const Geometry_> resol_;
  const Variables_ vars_;
  unsigned int rank_;
  boost::ptr_vector<Increment_> ensemblePerturbs_;
};

// ====================================================================================

template<typename MODEL>
Ensemble<MODEL>::Ensemble(const util::DateTime & validTime, const eckit::Configuration & conf)
  : config_(conf), validTime_(validTime), resol_(),
    vars_(eckit::LocalConfiguration(conf, "variables")),
    rank_(conf.getInt("members")),
    ensemblePerturbs_()
{
  Log::trace() << "Ensemble:contructor done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void Ensemble<MODEL>::linearize(const State_ & xb, const Geometry_ & resol) {
  ASSERT(xb.validTime() == validTime_);
  resol_.reset(new Geometry_(resol));
  State_ xblr(*resol_, xb);
  Accumulator<MODEL, Increment_, State_> bgmean(*resol_, vars_, validTime_);
  bgmean.accumul(1.0, xblr);
  const double rr = 1.0/static_cast<double>(rank_);

  std::vector<eckit::LocalConfiguration> confs;
  config_.get("state", confs);
  ASSERT(confs.size() == rank_);

  State_ xread(xblr);
  for (unsigned int jm = 0; jm < rank_; ++jm) {
    xread.read(confs[jm]);
    ASSERT(xread.validTime() == validTime_);

//  Ensemble will be centered around bg
    Increment_ * dx = new Increment_(*resol_, vars_, validTime_);
    dx->diff(xread, xblr);
    ensemblePerturbs_.push_back(dx);

//  Compute bg - ensemble mean
    bgmean.accumul(-rr, xread);
  }

// Re-center around mean instead of bg and rescale
  const double rk = 1.0 / sqrt((static_cast<double>(rank_) - 1.0));
  for (unsigned int jm = 0; jm < rank_; ++jm) {
    ensemblePerturbs_[jm] += bgmean;
    ensemblePerturbs_[jm] *= rk;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_ENSEMBLE_H_
