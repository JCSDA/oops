/*
 * (C) Copyright 2009-2016 ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_STATE_H_
#define OOPS_INTERFACE_STATE_H_

#include <string>

#include <boost/scoped_ptr.hpp>

#include "oops/base/PostProcessor.h"
#include "oops/interface/Geometry.h"
#include "oops/interface/Locations.h"
#include "oops/interface/ModelAtLocations.h"
#include "oops/interface/Variables.h"
#include "util/DateTime.h"
#include "util/ObjectCounter.h"
#include "util/Printable.h"
#include "util/Timer.h"

namespace eckit {
  class Configuration;
}

namespace oops {

/// Encapsulates the model state

// -----------------------------------------------------------------------------

template <typename MODEL>
class State : public util::Printable,
              private util::ObjectCounter<State<MODEL> > {
  typedef typename MODEL::State                 State_;
  typedef Geometry<MODEL>            Geometry_;
  typedef Locations<MODEL>           Locations_;
  typedef ModelAtLocations<MODEL>    ModelAtLocations_;
  typedef Variables<MODEL>           Variables_;

 public:
  static const std::string classname() {return "oops::State";}

/// Constructor, destructor
  State(const Geometry_ &, const eckit::Configuration &);
  State(const Geometry_ &, const State &);
  State(const State &);
  ~State();
  State & operator=(const State &);  // Is that used anywhere?

/// Interfacing
  State_ & state() {return *state_;}
  const State_ & state() const {return *state_;}

/// Interpolate to observation location
  void interpolate(const Locations_ &, ModelAtLocations_ &) const;

/// Time
  const util::DateTime validTime() const {return state_->validTime();}

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;
  Geometry_ geometry() const;

 protected:
/// Protected methods are for Accumulator. Could we find a better design?
  State(const Geometry_ &, const Variables_ &, const util::DateTime &);
  void zero();
  void accumul(const double &, const State &);

 private:
  void print(std::ostream &) const;
  boost::scoped_ptr<State_> state_;
};

// =============================================================================

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const Variables_ & vars,
                    const util::DateTime & time) : state_()
{
  Log::trace() << "State<MODEL>::State starting" << std::endl;
  util::Timer timer(classname(), "State");
  state_.reset(new State_(resol.geometry(), vars.variables(), time));
  Log::trace() << "State<MODEL>::State done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const eckit::Configuration & conf) : state_()
{
  Log::trace() << "State<MODEL>::State read starting" << std::endl;
  util::Timer timer(classname(), "State");
  state_.reset(new State_(resol.geometry(), conf));
  Log::trace() << "State<MODEL>::State read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const State & other) : state_()
{
  Log::trace() << "State<MODEL>::State interpolated starting" << std::endl;
  util::Timer timer(classname(), "State");
  state_.reset(new State_(resol.geometry(), *other.state_));
  Log::trace() << "State<MODEL>::State interpolated done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const State & other) : state_()
{
  Log::trace() << "State<MODEL>::State starting copy" << std::endl;
  util::Timer timer(classname(), "State");
  state_.reset(new State_(*other.state_));
  Log::trace() << "State<MODEL>::State copy done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::~State() {
  Log::trace() << "State<MODEL>::~State starting" << std::endl;
  util::Timer timer(classname(), "~State");
  state_.reset();
  Log::trace() << "State<MODEL>::~State done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL> & State<MODEL>::operator=(const State & rhs) {
  Log::trace() << "State<MODEL>::operator= starting" << std::endl;
  util::Timer timer(classname(), "operator=");
  *state_ = *rhs.state_;
  Log::trace() << "State<MODEL>::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::interpolate(const Locations_ & locs, ModelAtLocations_ & gom) const {
  Log::trace() << "State<MODEL>::interpolate starting" << std::endl;
  util::Timer timer(classname(), "interpolate");
  state_->interpolate(locs.locations(), gom.modelatlocations());
  Log::trace() << "State<MODEL>::interpolate done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::read(const eckit::Configuration & conf) {
  Log::trace() << "State<MODEL>::read starting" << std::endl;
  util::Timer timer(classname(), "read");
  state_->read(conf);
  Log::trace() << "State<MODEL>::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::write(const eckit::Configuration & conf) const {
  Log::trace() << "State<MODEL>::write starting" << std::endl;
  util::Timer timer(classname(), "write");
  state_->write(conf);
  Log::trace() << "State<MODEL>::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
double State<MODEL>::norm() const {
  Log::trace() << "State<MODEL>::norm starting" << std::endl;
  util::Timer timer(classname(), "norm");
  double zz = state_->norm();
  Log::trace() << "State<MODEL>::norm done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
Geometry<MODEL> State<MODEL>::geometry() const {
  Log::trace() << "State<MODEL>::geometry starting" << std::endl;
  util::Timer timer(classname(), "geometry");
  Geometry<MODEL> geom(state_->geometry());
  Log::trace() << "State<MODEL>::geometry done" << std::endl;
  return geom;
}

// -----------------------------------------------------------------------------


template<typename MODEL>
void State<MODEL>::print(std::ostream & os) const {
  Log::trace() << "State<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  os << *state_;
  Log::trace() << "State<MODEL>::print done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::zero() {
  Log::trace() << "State<MODEL>::zero starting" << std::endl;
  util::Timer timer(classname(), "zero");
  state_->zero();
  Log::trace() << "State<MODEL>::zero done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::accumul(const double & zz, const State & xx) {
  Log::trace() << "State<MODEL>::accumul starting" << std::endl;
  util::Timer timer(classname(), "accumul");
  state_->accumul(zz, *xx.state_);
  Log::trace() << "State<MODEL>::accumul done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

#endif  // OOPS_INTERFACE_STATE_H_
