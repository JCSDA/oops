/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_ASSIMILATION_STATE4D_H_
#define OOPS_ASSIMILATION_STATE4D_H_

#include <cmath>
#include <ostream>
#include <string>
#include <vector>

#include "eckit/config/LocalConfiguration.h"
#include "eckit/exception/Exceptions.h"

#include "oops/interface/Geometry.h"
#include "oops/interface/State.h"
#include "oops/util/Logger.h"
#include "oops/util/Printable.h"

namespace oops {

/// Four dimensional state
/*!
 *  The 4D state is mostly used as part of the VDA control variable.
 */

// -----------------------------------------------------------------------------
template<typename MODEL> class State4D : public util::Printable {
  typedef Geometry<MODEL>            Geometry_;
  typedef State<MODEL>               State_;

 public:
  static const std::string classname() {return "State4D";}

/// The arguments define the number of sub-windows and the resolution
  State4D(const Geometry_ &, const eckit::Configuration &);
  explicit State4D(const State_ &);

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;

/// Get model space control variable
  bool checkStatesNumber(const unsigned int nn) const {return state4d_.size() == nn;}
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
      if (config.has("member")) conf.set("member", config.getInt("member"));
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
State4D<MODEL>::State4D(const State_ & state3d) {
  state4d_.emplace_back(state3d);
  Log::trace() << "State4D constructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State4D<MODEL>::read(const eckit::Configuration & config) {
  // 4D state
  if (config.has("states")) {
    std::vector<eckit::LocalConfiguration> confs;
    config.get("states", confs);
    ASSERT(state4d_.size() == confs.size());
    for (size_t jj = 0; jj < state4d_.size(); ++jj) {
      state4d_[jj].read(confs[jj]);;
    }
  } else {
  // 3D state
    ASSERT(state4d_.size() == 1);
    state4d_[0].read(config);
  }
  for (size_t jj = 1; jj < state4d_.size(); ++jj) {
    ASSERT(state4d_[jj].variables() == state4d_[0].variables());
  }
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

template<typename MODEL>
double State4D<MODEL>::norm() const {
  double zn = 0.0;
  for (const State_ & state : state4d_) {
    double zz = state.norm();
    zn += zz * zz;
  }
  return sqrt(zn);
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_STATE4D_H_
