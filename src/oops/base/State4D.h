/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2020 UCAR
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_BASE_STATE4D_H_
#define OOPS_BASE_STATE4D_H_

#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"

#include "oops/base/Geometry.h"
#include "oops/interface/State.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

/// Four dimensional state (vector of 3D States)
template<typename MODEL> class State4D : public util::Printable {
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;

 public:
  static const std::string classname() {return "State4D";}

  /// The arguments define all states in 4D and their resolution
  State4D(const Geometry_ &, const eckit::Configuration &);

  /// I/O
  void write(const eckit::Configuration &) const;

  /// Get 3D model state
  size_t size() const {return state4d_.size();}
  State_ & operator[](const int ii) {return state4d_[ii];}
  const State_ & operator[](const int ii) const {return state4d_[ii];}

  Geometry_ geometry() const { return state4d_[0].geometry(); }
  const Variables & variables() const {return state4d_[0].variables();}
  const std::vector<util::DateTime> validTimes() const;

  /// Accumulator
  void zero();
  void accumul(const double &, const State4D &);

 private:
  void print(std::ostream &) const;

  std::vector<State_> state4d_;
};

// =============================================================================

template<typename MODEL>
State4D<MODEL>::State4D(const Geometry_ & resol, const eckit::Configuration & config) {
  Log::trace() << "State4D config : " << config << std::endl;
  // 4D state:
  if (config.has("states")) {
    std::vector<eckit::LocalConfiguration> confs;
    config.get("states", confs);
    state4d_.reserve(confs.size());
    for (auto & conf : confs) {
      state4d_.emplace_back(State_(resol, conf));
    }
  } else {
  // 3D state:
    state4d_.emplace_back(State_(resol, config));
  }
  for (size_t jj = 1; jj < state4d_.size(); ++jj) {
    ASSERT(state4d_[jj].variables() == state4d_[0].variables());
  }
  Log::trace() << "State4D constructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State4D<MODEL>::write(const eckit::Configuration & config) const {
  // 4D state
  if (config.has("states")) {
    std::vector<eckit::LocalConfiguration> confs;
    config.get("states", confs);
    ASSERT(state4d_.size() == confs.size());
    for (size_t jj = 0; jj < state4d_.size(); ++jj) {
      if (config.has("member")) confs[jj].set("member", config.getInt("member"));
      state4d_[jj].write(confs[jj]);
    }
  } else {
  // 3D state
    ASSERT(state4d_.size() == 1);
    state4d_[0].write(config);
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
const std::vector<util::DateTime> State4D<MODEL>::validTimes() const {
  std::vector<util::DateTime> times;
  times.reserve(state4d_.size());
  for (const State_ & state : state4d_) {
    times.push_back(state.validTime());
  }
  return times;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State4D<MODEL>::zero() {
  Log::trace() << "State4D<MODEL>::zero starting" << std::endl;
  for (State_ & state : state4d_) {
    state.zero();
  }
  Log::trace() << "State4D<MODEL>::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State4D<MODEL>::accumul(const double & zz, const State4D & xx) {
  Log::trace() << "State4D<MODEL>::accumul starting" << std::endl;
  ASSERT(xx.size() == state4d_.size());
  for (size_t jj = 0; jj < state4d_.size(); ++jj) {
    state4d_[jj].accumul(zz, xx[jj]);
  }
  Log::trace() << "State4D<MODEL>::accumul done" << std::endl;
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void State4D<MODEL>::print(std::ostream & outs) const {
  for (const State_ & state : state4d_) {
    outs << state << std::endl;
  }
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_BASE_STATE4D_H_
