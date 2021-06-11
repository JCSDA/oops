/*
 * (C) Copyright 2009-2016 ECMWF.
 * (C) Copyright 2017-2019 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation nor
 * does it submit to any jurisdiction.
 */

#ifndef OOPS_INTERFACE_STATE_H_
#define OOPS_INTERFACE_STATE_H_

#include <memory>
#include <string>
#include <vector>

#include "eckit/config/Configuration.h"
#include "oops/base/PostProcessor.h"
#include "oops/base/Variables.h"
#include "oops/interface/Geometry.h"
#include "oops/util/DateTime.h"
#include "oops/util/gatherPrint.h"
#include "oops/util/ObjectCounter.h"
#include "oops/util/Printable.h"
#include "oops/util/Serializable.h"
#include "oops/util/Timer.h"

namespace oops {

/// Encapsulates the model state

// -----------------------------------------------------------------------------

template <typename MODEL>
class State : public util::Printable,
              public util::Serializable,
              private util::ObjectCounter<State<MODEL> > {
  typedef typename MODEL::State      State_;
  typedef Geometry<MODEL>            Geometry_;

 public:
  static const std::string classname() {return "oops::State";}

/// Constructor, destructor
  State(const Geometry_ &, const Variables &, const util::DateTime &);
  State(const Geometry_ &, const eckit::Configuration &);
  State(const Geometry_ &, const State &);
  State(const State &);
  ~State();
  State & operator=(const State &);  // Is that used anywhere?

/// Interfacing
  State_ & state() {return *state_;}
  const State_ & state() const {return *state_;}

/// Time
  const util::DateTime validTime() const {return state_->validTime();}
  void updateTime(const util::Duration & dt) {state_->updateTime(dt);}

/// I/O and diagnostics
  void read(const eckit::Configuration &);
  void write(const eckit::Configuration &) const;
  double norm() const;  // Only for tests
  Geometry_ geometry() const;
  const Variables & variables() const;

/// Accumulator
  void zero();
  void accumul(const double &, const State &);

/// Serialize and deserialize
  size_t serialSize() const override;
  void serialize(std::vector<double> &) const override;
  void deserialize(const std::vector<double> &, size_t &) override;

 private:
  void print(std::ostream &) const override;
  std::unique_ptr<State_> state_;
  const eckit::mpi::Comm & commTime_;
};

// =============================================================================

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const Variables & vars,
                    const util::DateTime & time) : state_(), commTime_(resol.timeComm())
{
  Log::trace() << "State<MODEL>::State starting" << std::endl;
  util::Timer timer(classname(), "State");
  state_.reset(new State_(resol.geometry(), vars, time));
  this->setObjectSize(state_->serialSize()*sizeof(double));
  Log::trace() << "State<MODEL>::State done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const eckit::Configuration & conf)
  : state_(), commTime_(resol.timeComm())
{
  Log::trace() << "State<MODEL>::State read starting" << std::endl;
  util::Timer timer(classname(), "State");

  eckit::LocalConfiguration myconf;
  if (conf.has("states")) {
//  Parallel 4D state:
    std::vector<eckit::LocalConfiguration> confs;
    conf.get("states", confs);
    ASSERT(confs.size() == resol.timeComm().size());
    myconf = confs[resol.timeComm().rank()];
  } else {
//  3D state:
    myconf = eckit::LocalConfiguration(conf);
  }

  state_.reset(new State_(resol.geometry(), myconf));
  this->setObjectSize(state_->serialSize()*sizeof(double));
  Log::trace() << "State<MODEL>::State read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const Geometry_ & resol, const State & other)
  : state_(), commTime_(resol.timeComm())
{
  Log::trace() << "State<MODEL>::State interpolated starting" << std::endl;
  util::Timer timer(classname(), "State");
  state_.reset(new State_(resol.geometry(), *other.state_));
  this->setObjectSize(state_->serialSize()*sizeof(double));
  Log::trace() << "State<MODEL>::State interpolated done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
State<MODEL>::State(const State & other) : state_(), commTime_(other.commTime_)
{
  Log::trace() << "State<MODEL>::State starting copy" << std::endl;
  util::Timer timer(classname(), "State");
  state_.reset(new State_(*other.state_));
  this->setObjectSize(state_->serialSize()*sizeof(double));
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
  zz *= zz;
  commTime_.allReduceInPlace(zz, eckit::mpi::Operation::SUM);
  zz = sqrt(zz);
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
const Variables & State<MODEL>::variables() const {
  Log::trace() << "State<MODEL>::variables starting" << std::endl;
  util::Timer timer(classname(), "variables");
  return state_->variables();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
size_t State<MODEL>::serialSize() const {
  Log::trace() << "State<MODEL>::serialSize" << std::endl;
  util::Timer timer(classname(), "serialSize");
  return state_->serialSize();
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::serialize(std::vector<double> & vect) const {
  Log::trace() << "State<MODEL>::serialize starting" << std::endl;
  util::Timer timer(classname(), "serialize");
  state_->serialize(vect);
  Log::trace() << "State<MODEL>::serialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::deserialize(const std::vector<double> & vect, size_t & current) {
  Log::trace() << "State<MODEL>::State deserialize starting" << std::endl;
  util::Timer timer(classname(), "deserialize");
  state_->deserialize(vect, current);
  Log::trace() << "State<MODEL>::State deserialize done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL>
void State<MODEL>::print(std::ostream & os) const {
  Log::trace() << "State<MODEL>::print starting" << std::endl;
  util::Timer timer(classname(), "print");
  if (commTime_.size() > 1) {
    gatherPrint(os, *state_, commTime_);
  } else {
    os << *state_;
  }
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
