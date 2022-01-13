/*
 * (C) Copyright 2021-2021 UCAR.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 */

#pragma once

#include <memory>
#include <ostream>
#include <tuple>
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


// -----------------------------------------------------------------------------
/// Parameters for the State describing a coupled model state (used in ctor and read)
template<class... MODELs>
class StateCoupledParameters : public Parameters {
  OOPS_CONCRETE_PARAMETERS(StateCoupledParameters, Parameters)

  /// Type of tuples stored in the StateCoupledParameters
  using RequiredParametersTupleT =
        std::tuple<RequiredParameter<typename State<MODELs>::Parameters_>...>;
  /// Tuple that can be passed to the Parameter ctor
  using RequiredParameterInit = std::tuple<const char *, Parameters *>;

 public:
  /// Tuple of all State Parameters.
  RequiredParametersTupleT states{RequiredParameterInit(MODELs::name().c_str(), this) ... };
};

// -----------------------------------------------------------------------------
/// Parameters for the State describing how to output a coupled model state
template<class... MODELs>
class StateCoupledWriteParameters : public WriteParametersBase {
  OOPS_CONCRETE_PARAMETERS(StateCoupledWriteParameters, WriteParametersBase)

  /// Type of tuples stored in the StateCoupledWriteParameters
  using RequiredParametersTupleT =
        std::tuple<RequiredParameter<typename State<MODELs>::WriteParameters_>...>;
  /// Tuple that can be passed to the Parameter ctor
  using RequiredParameterInit = std::tuple<const char *, Parameters *>;

 public:
  /// Tuple of all Write state Parameters.
  RequiredParametersTupleT states{RequiredParameterInit(MODELs::name().c_str(), this) ... };
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
  State<MODEL1> & state1() {ASSERT(xx1_); return *xx1_;}
  State<MODEL2> & state2() {ASSERT(xx2_); return *xx2_;}
  const State<MODEL1> & state1() const {ASSERT(xx1_); return *xx1_;}
  const State<MODEL2> & state2() const {ASSERT(xx2_); return *xx2_;}

  /// Accessor to the state's valid time
  const util::DateTime validTime() const;
  /// Update state's time (increase by \p dt), used in the PseudoModel
  void updateTime(const util::Duration & dt);

  /// I/O and diagnostics
  void read(const Parameters_ &);
  void write(const WriteParameters_ &) const;
  double norm() const;
  std::shared_ptr<const GeometryCoupled_> geometry() const {return geom_;}
  const Variables & variables() const {return vars_;}

  // For accumulator
  void zero();
  void accumul(const double & zz, const StateCoupled & other);

  /// Serialize and deserialize
  size_t serialSize() const;
  void serialize(std::vector<double> & buf) const;
  void deserialize(const std::vector<double> & buf, size_t & ii);

 private:
  void print(std::ostream &) const;
  std::shared_ptr<const GeometryCoupled_> geom_;
  std::unique_ptr<State<MODEL1>> xx1_;
  std::unique_ptr<State<MODEL2>> xx2_;
  bool parallel_;
  Variables vars_;  // Add 'updateVars' method when implementing coupled variable change
};

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
StateCoupled<MODEL1, MODEL2>::StateCoupled(const GeometryCoupled_ & resol,
                                           const Variables & vars,
                                           const util::DateTime & time)
  : geom_(new GeometryCoupled_(resol)), xx1_(), xx2_(), parallel_(resol.isParallel()) {
  Log::trace() << "StateCoupled::StateCoupled starting" << std::endl;
  if (parallel_) {
    Log::debug() << "Parallel coupled states" << std::endl;
    if (resol.modelNumber() == 1) {
      xx1_ = std::make_unique<State<MODEL1>>(resol.geometry1(), vars, time);
      vars_ = xx1_->variables();
    }
    if (resol.modelNumber() == 2) {
      xx2_ = std::make_unique<State<MODEL2>>(resol.geometry2(), vars, time);
      vars_ = xx2_->variables();
    }
  } else {
    Log::debug() << "Sequential coupled states" << std::endl;
    xx1_ = std::make_unique<State<MODEL1>>(resol.geometry1(), vars, time);
    xx2_ = std::make_unique<State<MODEL2>>(resol.geometry2(), vars, time);
    vars_ = xx1_->variables();
    vars_ += xx2_->variables();
  }
  Log::trace() << "StateCoupled::StateCoupled done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
StateCoupled<MODEL1, MODEL2>::StateCoupled(const GeometryCoupled_ & resol,
                                           const Parameters_ & params)
  : geom_(new GeometryCoupled_(resol)), xx1_(), xx2_(), parallel_(resol.isParallel()) {
  Log::trace() << "StateCoupled::StateCoupled read starting" << std::endl;

  if (parallel_) {
    Log::debug() << "Parallel coupled states" << std::endl;
    if (resol.modelNumber() == 1) {
      xx1_ = std::make_unique<State<MODEL1>>(resol.geometry1(), std::get<0>(params.states));
      vars_ = xx1_->variables();
    }
    if (resol.modelNumber() == 2) {
      xx2_ = std::make_unique<State<MODEL2>>(resol.geometry2(), std::get<1>(params.states));
      vars_ = xx2_->variables();
    }
  } else {
    Log::debug() << "Sequential coupled states" << std::endl;
    xx1_ = std::make_unique<State<MODEL1>>(resol.geometry1(), std::get<0>(params.states));
    xx2_ = std::make_unique<State<MODEL2>>(resol.geometry2(), std::get<1>(params.states));
    vars_ = xx1_->variables();
    vars_ += xx2_->variables();
  }
  Log::trace() << "StateCoupled::StateCoupled read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
StateCoupled<MODEL1, MODEL2>::StateCoupled(const GeometryCoupled_ & resol,
                                           const StateCoupled & other)
  : geom_(new GeometryCoupled_(resol)), xx1_(), xx2_(), parallel_(resol.isParallel()),
    vars_(other.vars_) {
  Log::trace() << "StateCoupled::StateCoupled interpolated starting" << std::endl;
  ASSERT(parallel_ == other.parallel_);
  if (other.xx1_) xx1_ = std::make_unique<State<MODEL1>>(resol.geometry1(), *other.xx1_);
  if (other.xx2_) xx2_ = std::make_unique<State<MODEL2>>(resol.geometry2(), *other.xx2_);
  Log::trace() << "StateCoupled::StateCoupled interpolated done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
StateCoupled<MODEL1, MODEL2>::StateCoupled(const StateCoupled & other)
  : geom_(other.geom_), xx1_(), xx2_(), parallel_(other.parallel_), vars_(other.vars_) {
  Log::trace() << "StateCoupled::StateCoupled starting copy" << std::endl;
  if (other.xx1_) xx1_ = std::make_unique<State<MODEL1>>(*other.xx1_);
  if (other.xx2_) xx2_ = std::make_unique<State<MODEL2>>(*other.xx2_);
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
  if (xx1_) *xx1_ = *rhs.xx1_;
  if (xx2_) *xx2_ = *rhs.xx2_;
  Log::trace() << "StateCoupled::operator= done" << std::endl;
  return *this;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void StateCoupled<MODEL1, MODEL2>::read(const Parameters_ & params) {
  Log::trace() << "StateCoupled::read starting" << std::endl;
  if (xx1_) xx1_->read(std::get<0>(params.states));
  if (xx2_) xx2_->read(std::get<1>(params.states));
  Log::trace() << "StateCoupled::read done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void StateCoupled<MODEL1, MODEL2>::write(const WriteParameters_ & params) const {
  Log::trace() << "StateCoupled::write starting" << std::endl;
  if (xx1_) xx1_->write(std::get<0>(params.states));
  if (xx2_) xx2_->write(std::get<1>(params.states));
  Log::trace() << "StateCoupled::write done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
double StateCoupled<MODEL1, MODEL2>::norm() const {
  Log::trace() << "StateCoupled::norm starting" << std::endl;
  double zz = 0.0;
  if (parallel_) {
    if (xx1_) zz = xx1_->norm();
    if (xx2_) zz = xx2_->norm();
    geom_->getComm().allReduceInPlace(zz, eckit::mpi::Operation::SUM);
  } else {
    zz = xx1_->norm() + xx2_->norm();
  }
  return zz;
  Log::trace() << "StateCoupled::norm done" << std::endl;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
const util::DateTime StateCoupled<MODEL1, MODEL2>::validTime() const {
  util::DateTime dt;
  if (xx1_) dt = xx1_->validTime();
  if (xx2_) dt = xx2_->validTime();
  return dt;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void StateCoupled<MODEL1, MODEL2>::updateTime(const util::Duration & dt) {
  if (xx1_) xx1_->updateTime(dt);
  if (xx2_) xx2_->updateTime(dt);
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void StateCoupled<MODEL1, MODEL2>::zero() {
  if (xx1_) xx1_->zero();
  if (xx2_) xx2_->zero();
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void StateCoupled<MODEL1, MODEL2>::accumul(const double & zz, const StateCoupled & other) {
  if (xx1_) xx1_->accumul(zz, *other.xx1_);
  if (xx2_) xx2_->accumul(zz, *other.xx2_);
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
size_t StateCoupled<MODEL1, MODEL2>::serialSize() const {
  size_t ss = 0;
  if (xx1_) ss += xx1_->serialSize();
  if (xx2_) ss += xx2_->serialSize();
  return ss;
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void StateCoupled<MODEL1, MODEL2>::serialize(std::vector<double> & buf) const {
  if (xx1_) xx1_->serialize(buf);
  if (xx2_) xx2_->serialize(buf);
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void StateCoupled<MODEL1, MODEL2>::deserialize(const std::vector<double> & buf, size_t & ii) {
  if (xx1_) xx1_->deserialize(buf, ii);
  if (xx2_) xx2_->deserialize(buf, ii);
}

// -----------------------------------------------------------------------------

template<typename MODEL1, typename MODEL2>
void StateCoupled<MODEL1, MODEL2>::print(std::ostream & os) const {
  Log::trace() << "StateCoupled::print starting" << std::endl;
  if (parallel_ && geom_->localRank() == 0) {
    std::stringstream ss;
    if (xx1_) {
      ss << std::endl << "StateCoupled: " << MODEL1::name() << std::endl;
      ss << std::setprecision(os.precision()) << *xx1_ << std::endl;
    }
    if (xx2_) {
      ss << std::endl << "StateCoupled: " << MODEL2::name() << std::endl;
      ss << std::setprecision(os.precision()) << *xx2_ << std::endl;
    }
    util::gatherPrint(os, ss.str(), geom_->getComm());
  }

  if (!parallel_) {
    os << std::endl << "StateCoupled: " << MODEL1::name() << std::endl;
    os << *xx1_ << std::endl;
    os << std::endl << "StateCoupled: " << MODEL2::name() << std::endl;
    os << *xx2_;
  }
  Log::trace() << "StateCoupled::print done" << std::endl;
}
// -----------------------------------------------------------------------------

}  // namespace oops
