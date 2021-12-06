/*
 * (C) Copyright 2021-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <utility>
#include <vector>

#include "eckit/exception/Exceptions.h"

#include "oops/base/State.h"
#include "oops/base/Variables.h"
#include "oops/base/WriteParametersBase.h"
#include "oops/util/DateTime.h"
#include "oops/util/Duration.h"
#include "oops/util/Printable.h"

#include "oops/coupled/GeometryCoupled.h"

namespace oops {

/// Parameters for the State describing a coupled model state (used in ctor and read)
template <typename MODEL1, typename MODEL2>
class StateCoupledParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(StateCoupledParameters, Parameters)

  typedef typename State<MODEL1>::Parameters_ Parameters1_;
  typedef typename State<MODEL2>::Parameters_ Parameters2_;
 public:
  /// Parameters for State of MODEL1 and State of MODEL2
  RequiredParameter<Parameters1_> state1{MODEL1::name().c_str(), this};
  RequiredParameter<Parameters2_> state2{MODEL2::name().c_str(), this};
};

/// Parameters for the State describing how to output a coupled model state
template <typename MODEL1, typename MODEL2>
class StateCoupledWriteParameters : public WriteParametersBase {
  OOPS_CONCRETE_PARAMETERS(StateCoupledWriteParameters, WriteParametersBase)

  typedef typename State<MODEL1>::WriteParameters_ Parameters1_;
  typedef typename State<MODEL2>::WriteParameters_ Parameters2_;
 public:
  /// Parameters for State of MODEL1 and State of MODEL2
  RequiredParameter<Parameters1_> state1{MODEL1::name().c_str(), this};
  RequiredParameter<Parameters2_> state2{MODEL2::name().c_str(), this};
};

// -----------------------------------------------------------------------------
/// Coupled model state
template <typename MODEL1, typename MODEL2>
class StateCoupled : public util::Printable {
  typedef GeometryCoupled<MODEL1, MODEL2>  GeometryCoupled_;

 public:
  typedef StateCoupledParameters<MODEL1, MODEL2>      Parameters_;
  typedef StateCoupledWriteParameters<MODEL1, MODEL2> WriteParameters_;

  /// Constructor, destructor
  StateCoupled(const GeometryCoupled_ &, const Variables &, const util::DateTime &);
  StateCoupled(const GeometryCoupled_ &, const Parameters_ &);
  StateCoupled(const GeometryCoupled_ &, const StateCoupled &);
  StateCoupled(const StateCoupled &);
  virtual ~StateCoupled();
  StateCoupled & operator=(const StateCoupled &);

  /// Accessors to the coupled state components
  State<MODEL1> & state1() {return *xx1_;}
  State<MODEL2> & state2() {return *xx2_;}
  const State<MODEL1> & state1() const {return *xx1_;}
  const State<MODEL2> & state2() const {return *xx2_;}

  /// Accessor to the state's valid time
  const util::DateTime validTime() const
    {ASSERT(xx1_->validTime() == xx2_->validTime()); return xx1_->validTime();}
  /// Update state's time (increase by \p dt), used in the PseudoModel
  void updateTime(const util::Duration & dt) {xx1_->updateTime(dt); xx2_->updateTime(dt);}

  /// I/O and diagnostics
  void read(const Parameters_ &);
  void write(const WriteParameters_ &) const;
  double norm() const;
  std::shared_ptr<const GeometryCoupled_> geometry() const {return geom_;}
  const Variables & variables() const {return xx1_->variables();}

  // For accumulator
  void zero() {xx1_->zero(); xx2_->zero();}
  void accumul(const double & zz, const StateCoupled & other)
    {xx1_->accumul(zz, *other.xx1_); xx2_->accumul(zz, *other.xx2_);}

  /// Serialize and deserialize
  size_t serialSize() const {return xx1_->serialSize() + xx2_->serialSize();}
  void serialize(std::vector<double> & buf) const {xx1_->serialize(buf); xx2_->serialize(buf);}
  void deserialize(const std::vector<double> & buf, size_t & ii)
    {xx1_->deserialize(buf, ii); xx2_->deserialize(buf, ii);}

 private:
  void print(std::ostream &) const;
  std::shared_ptr<const GeometryCoupled_> geom_;
  std::unique_ptr<State<MODEL1>> xx1_;
  std::unique_ptr<State<MODEL2>> xx2_;
};

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
StateCoupled<MODEL1, MODEL2>::StateCoupled(const GeometryCoupled_ & resol,
                                           const Variables & vars,
                                           const util::DateTime & time)
  : geom_(new GeometryCoupled_(resol)), xx1_(), xx2_() {
  Log::trace() << "StateCoupled::StateCoupled starting" << std::endl;
  xx1_ = std::make_unique<State<MODEL1>>(resol.geometry1(), vars, time);
  xx2_ = std::make_unique<State<MODEL2>>(resol.geometry2(), vars, time);
  Log::trace() << "StateCoupled::StateCoupled done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
StateCoupled<MODEL1, MODEL2>::StateCoupled(const GeometryCoupled_ & resol,
                                           const Parameters_ & params)
  : geom_(new GeometryCoupled_(resol)), xx1_(), xx2_() {
  Log::trace() << "StateCoupled::StateCoupled read starting" << std::endl;

  xx1_ = std::make_unique<State<MODEL1>>(resol.geometry1(), params.state1);
  xx2_ = std::make_unique<State<MODEL2>>(resol.geometry2(), params.state2);

  Log::trace() << "StateCoupled::StateCoupled read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
StateCoupled<MODEL1, MODEL2>::StateCoupled(const GeometryCoupled_ & resol,
                                           const StateCoupled & other)
  : geom_(new GeometryCoupled_(resol)), xx1_(), xx2_()
{
  Log::trace() << "StateCoupled::StateCoupled interpolated starting" << std::endl;
  xx1_ = std::make_unique<State<MODEL1>>(resol.geometry1(), *other.xx1_);
  xx2_ = std::make_unique<State<MODEL2>>(resol.geometry2(), *other.xx2_);
  Log::trace() << "StateCoupled::StateCoupled interpolated done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
StateCoupled<MODEL1, MODEL2>::StateCoupled(const StateCoupled & other)
  : geom_(other.geom_), xx1_(), xx2_() {
  Log::trace() << "StateCoupled::StateCoupled starting copy" << std::endl;
  xx1_ = std::make_unique<State<MODEL1>>(*other.xx1_);
  xx2_ = std::make_unique<State<MODEL2>>(*other.xx2_);
  Log::trace() << "StateCoupled::StateCoupled copy done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
StateCoupled<MODEL1, MODEL2>::~StateCoupled() {
  Log::trace() << "StateCoupled::~StateCoupled starting" << std::endl;
  xx1_.reset();
  xx2_.reset();
  Log::trace() << "StateCoupled::~StateCoupled done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
StateCoupled<MODEL1, MODEL2> & StateCoupled<MODEL1, MODEL2>::operator=(const StateCoupled & rhs) {
  Log::trace() << "StateCoupled::operator= starting" << std::endl;
  *xx1_ = *rhs.xx1_;
  *xx2_ = *rhs.xx2_;
  Log::trace() << "StateCoupled::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void StateCoupled<MODEL1, MODEL2>::read(const Parameters_ & params) {
  Log::trace() << "StateCoupled::read starting" << std::endl;
  xx1_->read(params.state1);
  xx2_->read(params.state2);
  Log::trace() << "StateCoupled::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void StateCoupled<MODEL1, MODEL2>::write(const WriteParameters_ & params) const {
  Log::trace() << "StateCoupled::write starting" << std::endl;
  xx1_->write(params.state1);
  xx2_->write(params.state2);
  Log::trace() << "StateCoupled::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
double StateCoupled<MODEL1, MODEL2>::norm() const {
  Log::trace() << "StateCoupled::norm starting" << std::endl;
  const double zz = xx1_->norm() + xx2_->norm();
  Log::trace() << "StateCoupled::norm done" << std::endl;
  return zz;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void StateCoupled<MODEL1, MODEL2>::print(std::ostream & os) const {
  Log::trace() << "StateCoupled::print starting" << std::endl;
  os << std::endl << "StateCoupled: " << MODEL1::name();
  os << *xx1_;
  os << std::endl << "StateCoupled: " << MODEL2::name();
  os << *xx2_;
  Log::trace() << "StateCoupled::print done" << std::endl;
}

// -----------------------------------------------------------------------------

}  // namespace oops

