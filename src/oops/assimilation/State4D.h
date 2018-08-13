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
#include <boost/foreach.hpp>
#include <boost/ptr_container/ptr_vector.hpp>

#include "eckit/config/LocalConfiguration.h"
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
  State4D(const eckit::Configuration &, const Geometry_ &);
  explicit State4D(const State4D &);
  ~State4D();

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;

/// Get model space control variable
  bool checkStatesNumber(const unsigned int nn) const {return state4d_.size() == nn;}
  size_t size() const {return state4d_.size();}
  State_ & operator[](const int ii) {return state4d_[ii];}
  const State_ & operator[](const int ii) const {return state4d_[ii];}

 private:
  State4D & operator= (const State4D &);  // No assignment
  void print(std::ostream &) const;

  boost::ptr_vector<State_> state4d_;
};

// =============================================================================

template<typename MODEL>
State4D<MODEL>::State4D(const eckit::Configuration & config, const Geometry_ & resol) {
  std::vector<eckit::LocalConfiguration> files;
  config.get("state", files);
  Log::debug() << "State4D: reading " << files.size() << " states." << std::endl;

  for (size_t jsub = 0; jsub < files.size(); ++jsub) {
    Log::debug() << "State4D:reading" << files[jsub] << std::endl;
    State_ * js = new State_(resol, files[jsub]);
    Log::debug() << "State4D:State4D: read bg at " << js->validTime() << std::endl;
    state4d_.push_back(js);
  }
  Log::trace() << "State4D constructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State4D<MODEL>::State4D(const State4D & other) {
  BOOST_FOREACH(const State_ & js, other.state4d_)
    state4d_.push_back(new State_(js));
  Log::trace() << "State4D copied." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State4D<MODEL>::~State4D() {
  Log::trace() << "State4D destructed." << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State4D<MODEL>::read(const eckit::Configuration & config) {
  std::vector<eckit::LocalConfiguration> confs;
  config.get("state", confs);
  ASSERT(state4d_.size() == confs.size());
  unsigned jsub = 0;
  BOOST_FOREACH(State_ & js, state4d_) {
    Log::debug() << "State4D:read" << confs[jsub] << std::endl;
    js.read(confs[jsub]);
    ++jsub;
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State4D<MODEL>::write(const eckit::Configuration & config) const {
  std::vector<eckit::LocalConfiguration> confs;
  config.get("state", confs);
  ASSERT(state4d_.size() == confs.size());
  unsigned int jsub = 0;
  BOOST_FOREACH(const State_ & js, state4d_) {
    Log::debug() << "State4D:write" << confs[jsub] << std::endl;
    js.write(confs[jsub]);
    ++jsub;
  }
}

// -----------------------------------------------------------------------------

template <typename MODEL>
void State4D<MODEL>::print(std::ostream & outs) const {
  BOOST_FOREACH(const State_ & js, state4d_) {
    outs << js << std::endl;
  }
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double State4D<MODEL>::norm() const {
  double zn = 0.0;
  BOOST_FOREACH(const State_ & js, state4d_) {
    double zz = js.norm();
    zn += zz * zz;
  }
  return sqrt(zn);
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_ASSIMILATION_STATE4D_H_
